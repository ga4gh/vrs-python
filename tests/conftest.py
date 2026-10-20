import os
from collections.abc import Iterator
from pathlib import Path

import hgvs.dataproviders.uta
import pytest
from biocommons.seqrepo import SeqRepo
from hgvs_cache_json import json_to_pickle
from offline_guard import block_network, require_recorded_data

from ga4gh.vrs.dataproxy import SeqRepoDataProxy, SeqRepoRESTDataProxy

HGVS_CACHE_MODES = ("learn", "run")
HGVS_CACHE_JSON = Path(__file__).parent / "data" / "hgvs_cache.json"
# Run mode never connects to UTA. This URL refuses connections, and being on localhost
# it also keeps hgvs from logging that tests are using remote data.
HGVS_OFFLINE_UTA_URL = "postgresql://anonymous@localhost:1/uta/uta_offline"


def _hgvs_cache_mode() -> str:
    """Return `VRS_HGVS_CACHE_MODE` (default `run`)"""
    return os.environ.get("VRS_HGVS_CACHE_MODE", "run")


def pytest_configure(config: pytest.Config) -> None:
    """Check the hgvs cache settings once, before any test runs"""
    mode = _hgvs_cache_mode()
    if mode not in HGVS_CACHE_MODES:
        msg = f"VRS_HGVS_CACHE_MODE must be one of {HGVS_CACHE_MODES}, got {mode!r}"
        raise pytest.UsageError(msg)
    if mode == "learn" and config.getoption("--vcr-record") == "none":
        msg = "HGVS recording requires live connections; use --disable-vcr without --vcr-record=none."
        raise pytest.UsageError(msg)
    if mode == "learn" and not os.environ.get("VRS_HGVS_CACHE_FILE"):
        msg = (
            "VRS_HGVS_CACHE_FILE must name the pickle file to record into in learn mode"
        )
        raise pytest.UsageError(msg)


@pytest.fixture(scope="session")
def hgvs_cache_mode() -> str:
    """Return the hgvs cache mode tests run in: `run` or `learn`"""
    return _hgvs_cache_mode()


@pytest.fixture(scope="session")
def network_blocked(request: pytest.FixtureRequest) -> Iterator[bool]:
    """Block live connections before provider setup when replaying VCR cassettes"""
    blocked = request.config.getoption("--vcr-record") == "none"
    with pytest.MonkeyPatch.context() as mp:
        if blocked:
            block_network(mp)
        yield blocked


@pytest.fixture(scope="session", autouse=True)
def hgvs_cached_data_provider(
    hgvs_cache_mode: str,
    network_blocked: bool,  # noqa: ARG001  # Install blockers before provider setup.
    tmp_path_factory: pytest.TempPathFactory,
) -> Iterator[hgvs.dataproviders.uta.UTABase]:
    """Share one provider for HGVS replay or recording.

    Replay (`run`, the default) loads the committed JSON via a temporary pickle and
    fails immediately on missing data. `make record-hgvs-cache` uses `learn` with
    VRS_HGVS_CACHE_FILE and live services to record a replacement.
    """
    if hgvs_cache_mode == "learn":
        cache_file = os.environ["VRS_HGVS_CACHE_FILE"]
    else:
        cache_file = str(tmp_path_factory.mktemp("hgvs_cache") / "hgvs_cache.pkl")
        json_to_pickle(HGVS_CACHE_JSON, Path(cache_file))

    db_url = HGVS_OFFLINE_UTA_URL if hgvs_cache_mode == "run" else None
    provider = hgvs.dataproviders.uta.connect(
        db_url=db_url, mode=hgvs_cache_mode, cache=cache_file
    )
    try:
        with pytest.MonkeyPatch.context() as mp:
            if hgvs_cache_mode == "run":
                # In run mode, these wrappers raise HGVSDataNotAvailableError only
                # for cache misses; ordinary domain errors occur above this layer.
                for name, method in vars(provider).items():
                    if hasattr(method, "cache_info"):
                        mp.setattr(provider, name, require_recorded_data(method))
            mp.setattr(
                hgvs.dataproviders.uta, "connect", lambda *_args, **_kwargs: provider
            )
            # HgvsTools.close() closes its data provider, but this one is shared.
            mp.setattr(provider, "close", lambda: None)
            yield provider
    finally:
        provider.close()


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
