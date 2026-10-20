import os
from collections.abc import Iterator
from pathlib import Path

import hgvs.dataproviders.uta
import pytest
from biocommons.seqrepo import SeqRepo
from hgvs_cache_json import json_to_pickle
from offline_guard import OfflineGuard

from ga4gh.vrs.dataproxy import SeqRepoDataProxy, SeqRepoRESTDataProxy

HGVS_CACHE_MODES = ("learn", "run", "verify")
HGVS_CACHE_JSON = Path(__file__).parent / "data" / "hgvs_cache.json"
# Run mode never connects to UTA. This URL refuses connections, and being on localhost
# it also keeps hgvs from logging that tests are using remote data.
HGVS_OFFLINE_UTA_URL = "postgresql://anonymous@localhost:1/uta/uta_offline"


def _hgvs_cache_mode() -> str | None:
    """Return `VRS_HGVS_CACHE_MODE` (default `run`), or None when it's empty"""
    return os.environ.get("VRS_HGVS_CACHE_MODE", "run") or None


def pytest_configure() -> None:
    """Check the hgvs cache settings once, before any test runs"""
    mode = _hgvs_cache_mode()
    if mode is not None and mode not in HGVS_CACHE_MODES:
        msg = f"VRS_HGVS_CACHE_MODE must be one of {HGVS_CACHE_MODES} or empty, got {mode!r}"
        raise pytest.UsageError(msg)
    if mode == "learn" and not os.environ.get("VRS_HGVS_CACHE_FILE"):
        msg = (
            "VRS_HGVS_CACHE_FILE must name the pickle file to record into in learn mode"
        )
        raise pytest.UsageError(msg)


@pytest.fixture(scope="session")
def hgvs_cache_mode() -> str | None:
    """Return the hgvs cache mode tests run in: `run`, `learn`, `verify`, or None"""
    return _hgvs_cache_mode()


@pytest.fixture(scope="session", autouse=True)
def hgvs_cached_data_provider(
    hgvs_cache_mode: str | None,
    tmp_path_factory: pytest.TempPathFactory,
) -> Iterator[hgvs.dataproviders.uta.UTABase | None]:
    """Serve hgvs data provider lookups (UTA queries and hgvs sequence fetches) from a
    recorded cache, so tests don't need a UTA database or network access.

    Every `HgvsTools` gets this one data provider, so in learn mode they cannot
    overwrite each other's cache entries. `VRS_HGVS_CACHE_MODE` is `run` (default),
    `learn`, or `verify`, or empty to query UTA directly without the cache. In run
    mode, `offline_check` fails any test during which hgvs misses the cache.

    hgvs keeps its cache as a pickle, so in run and verify modes the committed
    `tests/data/hgvs_cache.json` is converted to a temporary pickle for hgvs to load
    (see `hgvs_cache_json.py`). Learn mode records a new cache, starting empty, into
    the pickle file named by `VRS_HGVS_CACHE_FILE`. To re-record, run
    `make record-hgvs-cache` with UTA_DB_URL pointing at a UTA instance and a
    seqrepo-rest-service running at SEQREPO_REST_URL; it converts the recorded pickle
    to `tests/data/hgvs_cache.json`.

    Learn and verify modes, and an empty mode, query UTA and fetch sequences live. The
    VCR cassettes don't contain those hgvs sequence requests (to HGVS_SEQREPO_URL), so
    run these modes with `pytest --disable-vcr`.
    """
    if hgvs_cache_mode is None:
        yield None
        return
    if hgvs_cache_mode == "learn":
        cache_file = os.environ["VRS_HGVS_CACHE_FILE"]
    else:
        cache_file = str(tmp_path_factory.mktemp("hgvs_cache") / "hgvs_cache.pkl")
        json_to_pickle(HGVS_CACHE_JSON, Path(cache_file))

    db_url = HGVS_OFFLINE_UTA_URL if hgvs_cache_mode == "run" else None
    provider = hgvs.dataproviders.uta.connect(
        db_url=db_url, mode=hgvs_cache_mode, cache=cache_file
    )
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(
            hgvs.dataproviders.uta, "connect", lambda *_args, **_kwargs: provider
        )
        # HgvsTools.close() closes its data provider, but this one is shared
        mp.setattr(provider, "close", lambda: None)
        yield provider
    provider.close()


@pytest.fixture(scope="session")
def offline_guard(
    request: pytest.FixtureRequest,
    hgvs_cache_mode: str | None,
    hgvs_cached_data_provider: hgvs.dataproviders.uta.UTABase | None,
) -> Iterator[OfflineGuard]:
    """Guard the session against going past recorded data (see `offline_guard.py`)

    With `--vcr-record=none`, network connections are blocked; other record modes
    record cassettes, which needs the network. With the hgvs cache in run mode, hgvs
    cache misses are recorded.
    """
    guard = OfflineGuard()
    with pytest.MonkeyPatch.context() as mp:
        if request.config.getoption("--vcr-record") == "none":
            guard.block_network(mp)
        if hgvs_cache_mode == "run":
            guard.record_hgvs_cache_misses(mp, hgvs_cached_data_provider)
        yield guard


@pytest.fixture(autouse=True)
def offline_check(offline_guard: OfflineGuard) -> Iterator[None]:
    """Fail a test, at teardown, if a network connection was opened or the hgvs cache
    was missed since the previous test's check, which includes the setup and teardown
    of wider-scoped fixtures
    """
    yield
    if violations := offline_guard.violations[:]:
        offline_guard.violations.clear()
        pytest.fail(
            "This test, or a fixture set up or torn down before it, went past its "
            "recorded data. Record its HTTP requests in a VCR cassette "
            "(`@pytest.mark.vcr`), and re-record hgvs data with "
            "`make record-hgvs-cache`.\n" + "\n".join(violations)
        )


def remove_request_headers(request):
    """Remove all headers from VCR request before recording."""
    request.headers = {}
    return request


def remove_response_headers(response):
    """Remove all headers from VCR response before recording."""
    response["headers"] = {}
    return response


@pytest.fixture(scope="module")
def vcr_config():
    """Configure VCR to filter out headers from cassettes."""
    return {
        "before_record_request": remove_request_headers,
        "before_record_response": remove_response_headers,
        "decode_compressed_response": True,
    }


@pytest.fixture(scope="session")
def dataproxy():
    sr = SeqRepo(
        root_dir=os.environ.get("SEQREPO_ROOT_DIR", "/usr/local/share/seqrepo/latest")
    )
    return SeqRepoDataProxy(sr)


@pytest.fixture(scope="session")
def rest_dataproxy():
    return SeqRepoRESTDataProxy(
        base_url=os.environ.get("SEQREPO_REST_URL", "http://localhost:5000/seqrepo"),
        disable_healthcheck=True,
    )


# See https://github.com/ga4gh/vrs-python/issues/24
# @pytest.fixture(autouse=True)
# def setup_doctest(doctest_namespace, tlr):
#     doctest_namespace["tlr"] = tlr
