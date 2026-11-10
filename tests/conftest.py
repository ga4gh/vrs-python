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
        help="Record HGVS data to a pickle file (requires --disable-vcr).",
    )


def pytest_configure(config: pytest.Config) -> None:
    """Reject recording settings that would block or replay live requests"""
    if config.getoption("--record-hgvs-cache") and not config.getoption(
        "--disable-vcr"
    ):
        msg = "HGVS recording requires --disable-vcr."
        raise pytest.UsageError(msg)


def _forbidden_connection(_sock: socket.socket, address: object) -> NoReturn:
    pytest.fail(
        f"Live connection to {address!r} is disabled. "
        "Record HTTP requests in a VCR cassette (@pytest.mark.vcr) with --vcr-record=once."
    )


def _forbidden_postgres(*_args: object, **_kwargs: object) -> NoReturn:
    pytest.fail("Live PostgreSQL connections are disabled. Run make record-hgvs-cache.")


def _require_recorded_data(method: Callable[..., Any]) -> Callable[..., Any]:
    """Make cache misses fatal even inside production HGVS fallback handlers"""

    @functools.wraps(method)
    def wrapped(*args: object, **kwargs: object) -> Any:  # noqa: ANN401
        try:
            return method(*args, **kwargs)
        except HGVSDataNotAvailableError as exc:
            pytest.fail(f"Missing HGVS recording: {exc}\nRun make record-hgvs-cache.")

    return wrapped


@pytest.fixture(scope="session", autouse=True)
def hgvs_cached_data_provider(
    request: pytest.FixtureRequest,
    tmp_path_factory: pytest.TempPathFactory,
) -> Iterator[hgvs.dataproviders.uta.UTABase]:
    """Share one provider so recording instances cannot overwrite each other's entries.

    Replay converts committed JSON to a temporary pickle and blocks live connections
    unless HTTP recording or --disable-vcr explicitly enables them.
    """
    recording = request.config.getoption("--record-hgvs-cache")
    cache_file = recording or tmp_path_factory.mktemp("hgvs_cache") / "hgvs_cache.pkl"
    if not recording:
        json_to_pickle(HGVS_CACHE_JSON, cache_file)

    with pytest.MonkeyPatch.context() as mp:
        if request.config.getoption(
            "--vcr-record"
        ) == "none" and not request.config.getoption("--disable-vcr"):
            # VCR replays without connecting. libpq bypasses Python sockets, so
            # psycopg2 needs its own patch. pytest.fail escapes except Exception.
            mp.setattr(socket.socket, "connect", _forbidden_connection)
            mp.setattr(socket.socket, "connect_ex", _forbidden_connection)
            mp.setattr(psycopg2, "connect", _forbidden_postgres)
        provider = hgvs.dataproviders.uta.connect(
            db_url=None if recording else HGVS_OFFLINE_UTA_URL,
            mode="learn" if recording else "run",
            cache=str(cache_file),
        )
        close = provider.close
        try:
            if not recording:
                # In run mode, these methods raise HGVSDataNotAvailableError only
                # for cache misses; ordinary domain errors occur above this layer.
                for name, method in vars(provider).items():
                    if hasattr(method, "cache_info"):
                        mp.setattr(provider, name, _require_recorded_data(method))
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
