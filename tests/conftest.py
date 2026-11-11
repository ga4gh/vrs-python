import functools
import os
import socket
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any, NoReturn

import hgvs.dataproviders.uta
import psycopg2
import pytest
from biocommons.seqrepo import SeqRepo
from hgvs.decorators.lru_cache import _make_key
from hgvs.exceptions import HGVSDataNotAvailableError
from hgvs_cache_json import json_to_pickle

from ga4gh.vrs.dataproxy import SeqRepoDataProxy, SeqRepoRESTDataProxy

HGVS_CACHE_JSON = Path(__file__).parent / "data" / "hgvs_cache.json"
# Run mode never connects to UTA. This URL refuses connections, and being on localhost
# it also keeps hgvs from logging that tests are using remote data.
HGVS_OFFLINE_UTA_URL = "postgresql://anonymous@localhost:1/uta/uta_offline"


def pytest_addoption(parser: pytest.Parser) -> None:
    """Expose explicit HGVS recording; ordinary test runs replay committed data"""
    parser.addoption(
        "--record-hgvs-cache",
        metavar="PATH",
        type=Path,
        help=(
            "Record live hgvs lookups into the pickle file PATH, adding to it if it "
            "exists (requires --disable-vcr)."
        ),
    )


def pytest_configure(config: pytest.Config) -> None:
    """Reject --record-hgvs-cache settings that would fail only after tests start"""
    record_path: Path | None = config.getoption("--record-hgvs-cache")
    if record_path is None:
        return
    if not config.getoption("--disable-vcr"):
        msg = "HGVS recording requires --disable-vcr."
        raise pytest.UsageError(msg)
    # hgvs opens the file only after the first cache miss, so check it upfront.
    if record_path.exists() and not record_path.is_file():
        msg = f"--record-hgvs-cache must name a regular file: {record_path}"
        raise pytest.UsageError(msg)
    if not record_path.parent.is_dir():
        msg = f"--record-hgvs-cache directory does not exist: {record_path.parent}"
        raise pytest.UsageError(msg)
    # Without it, hgvs silently records from its built-in public UTA.
    if not os.environ.get("UTA_DB_URL"):
        msg = "HGVS recording requires UTA_DB_URL to name the UTA instance to record from."
        raise pytest.UsageError(msg)


def _forbidden_connection(_sock: socket.socket, address: object) -> NoReturn:
    pytest.fail(
        f"Live connection to {address!r} is disabled. "
        "Record HTTP requests in a VCR cassette (@pytest.mark.vcr) with --vcr-record=once."
    )


def _forbidden_postgres(*_args: object, **_kwargs: object) -> NoReturn:
    pytest.fail("Live PostgreSQL connections are disabled. Run make record-hgvs-cache.")


@pytest.fixture(scope="session", autouse=True)
def block_live_connections(request: pytest.FixtureRequest) -> Iterator[bool]:
    """Fail any test that opens a live connection while VCR replays cassettes.

    Active with `--vcr-record=none` (the default) unless `--disable-vcr` is given.
    Yields whether connections are blocked.
    """
    blocked = request.config.getoption(
        "--vcr-record"
    ) == "none" and not request.config.getoption("--disable-vcr")
    with pytest.MonkeyPatch.context() as mp:
        if blocked:
            # VCR replays without connecting. libpq bypasses Python sockets, so
            # psycopg2 needs its own patch. pytest.fail escapes except Exception.
            mp.setattr(socket.socket, "connect", _forbidden_connection)
            mp.setattr(socket.socket, "connect_ex", _forbidden_connection)
            mp.setattr(psycopg2, "connect", _forbidden_postgres)
        yield blocked


def _wrap_cached_method(
    provider: hgvs.dataproviders.uta.UTABase,
    method: Callable[..., Any],
    recording: bool,
) -> Callable[..., Any]:
    """Record and replay "no data" results, and fail on lookups never recorded.

    hgvs caches only returned values. Without this, a lookup that raised
    HGVSDataNotAvailableError would be missing from the recording, and replaying it
    would look the same as a cache miss. A recorded result is stored as the exception
    itself, which hgvs's cache returns like any other value.

    hgvs reports a SeqRepo outage as HGVSDataNotAvailableError too. An outage while
    recording fails tests, so no JSON is written, and any that slips through shows in
    the JSON diff.
    """
    name: str = method.__wrapped__.__name__  # hgvs keys cache entries by this name

    @functools.wraps(method)
    def wrapped(*args: object, **kwargs: object) -> Any:
        try:
            result = method(*args, **kwargs)
        except HGVSDataNotAvailableError as exc:
            if not recording:
                # Replay never queries a data source, so only a cache miss raises.
                # pytest.fail escapes the production handlers that catch this error.
                pytest.fail(
                    f"Missing HGVS recording: {exc}\nRun make record-hgvs-cache."
                )
            # Keep the message: hgvs's Normalizer.is_valid_pos reads it.
            provider.cache[_make_key(name, args, kwargs, False, ())] = (
                HGVSDataNotAvailableError(str(exc))
            )
            provider.cache.sync()
            raise
        if isinstance(result, HGVSDataNotAvailableError):
            # A fresh instance, so tracebacks don't pile up on the cached one.
            raise HGVSDataNotAvailableError(str(result))
        return result

    return wrapped


@pytest.fixture(scope="session", autouse=True)
def hgvs_cached_data_provider(
    request: pytest.FixtureRequest,
    tmp_path_factory: pytest.TempPathFactory,
    block_live_connections: bool,  # noqa: ARG001  # Block before provider setup.
) -> Iterator[hgvs.dataproviders.uta.UTABase]:
    """Share one provider so recording instances cannot overwrite each other's entries.

    Replay converts committed JSON to a temporary pickle; recording
    (`--record-hgvs-cache`) writes live data to the given pickle.
    """
    record_path: Path | None = request.config.getoption("--record-hgvs-cache")
    recording = record_path is not None
    cache_file = record_path or tmp_path_factory.mktemp("hgvs_cache") / "hgvs_cache.pkl"
    if not recording:
        json_to_pickle(HGVS_CACHE_JSON, cache_file)

    with pytest.MonkeyPatch.context() as mp:
        provider = hgvs.dataproviders.uta.connect(
            db_url=os.environ["UTA_DB_URL"] if recording else HGVS_OFFLINE_UTA_URL,
            mode="learn" if recording else "run",
            cache=str(cache_file),
        )
        close = provider.close
        try:
            for name, method in vars(provider).items():
                if hasattr(method, "cache_info"):
                    mp.setattr(
                        provider, name, _wrap_cached_method(provider, method, recording)
                    )
            mp.setattr(hgvs.dataproviders.uta, "connect", lambda *_a, **_kw: provider)
            # HgvsTools.close() must not close the shared provider between tests.
            mp.setattr(provider, "close", lambda: None)
            yield provider
        finally:
            close()


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
