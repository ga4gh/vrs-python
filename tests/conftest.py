import os
import pickle
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import hgvs.dataproviders.uta
import pytest
from biocommons.seqrepo import SeqRepo
from hgvs.decorators.lru_cache import _make_key

from ga4gh.vrs.dataproxy import SeqRepoDataProxy, SeqRepoRESTDataProxy

HGVS_CACHE_MODES = ("learn", "run", "verify")
HGVS_CACHE_FILE = Path(__file__).parent / "data" / "hgvs_cache.pkl"
# The only classes a recorded hgvs cache file is made of: its keys (_HashedSeq), UTA
# rows (DictRow, with an OrderedDict column index), and what pickle rebuilds them with
HGVS_CACHE_ALLOWED_CLASSES = {
    ("builtins", "list"),
    ("collections", "OrderedDict"),
    ("copyreg", "_reconstructor"),
    ("hgvs.decorators.lru_cache", "_HashedSeq"),
    ("psycopg2.extras", "DictRow"),
}


class HgvsCacheUnpickler(pickle.Unpickler):
    """Unpickler that refuses any class not in `HGVS_CACHE_ALLOWED_CLASSES`"""

    def find_class(self, module: str, name: str) -> Any:  # noqa: ANN401
        if (module, name) not in HGVS_CACHE_ALLOWED_CLASSES:
            msg = f"hgvs cache file references a disallowed class: {module}.{name}"
            raise pickle.UnpicklingError(msg)
        return super().find_class(module, name)


def check_hgvs_cache_file(cache_file: str) -> None:
    """Check an hgvs cache file before hgvs loads it with an unrestricted `pickle.load`

    The file may only reference the classes in `HGVS_CACHE_ALLOWED_CLASSES`, and its
    keys must be in the format this hgvs version builds. Every data provider looks up
    `schema_version` when it's created, so a recorded cache always has that key.
    """
    with Path(cache_file).open("rb") as f:
        cache: dict = HgvsCacheUnpickler(f).load()  # cache keys -> recorded results
    if _make_key("schema_version", (), {}, False, ()) not in cache:
        msg = (
            f"{cache_file} has no schema_version entry in this hgvs version's cache key "
            "format; re-record it with `make record-hgvs-cache`"
        )
        raise pytest.UsageError(msg)


@pytest.fixture(scope="session", autouse=True)
def hgvs_cached_data_provider() -> Iterator[hgvs.dataproviders.uta.UTABase | None]:
    """Serve hgvs data provider lookups (UTA queries and hgvs sequence fetches) from a
    recorded cache, so tests don't need a UTA database or network access.

    Every `HgvsTools` gets this one data provider, so in learn mode they cannot
    overwrite each other's cache entries. `VRS_HGVS_CACHE_MODE` is `run` (default),
    `learn`, or `verify`, or empty to query UTA directly without the cache, and
    `VRS_HGVS_CACHE_FILE` overrides the cache file. In run mode, a lookup missing from
    the cache raises `HGVSDataNotAvailableError`. To re-record, run
    `make record-hgvs-cache` with UTA_DB_URL pointing at a UTA instance and a
    seqrepo-rest-service running at SEQREPO_REST_URL.
    """
    mode = os.environ.get("VRS_HGVS_CACHE_MODE", "run")
    if not mode:
        yield None
        return
    if mode not in HGVS_CACHE_MODES:
        msg = f"VRS_HGVS_CACHE_MODE must be one of {HGVS_CACHE_MODES} or empty, got {mode!r}"
        raise pytest.UsageError(msg)
    cache_file = os.environ.get("VRS_HGVS_CACHE_FILE") or str(HGVS_CACHE_FILE)
    if Path(cache_file).exists():
        check_hgvs_cache_file(cache_file)

    provider = hgvs.dataproviders.uta.connect(mode=mode, cache=cache_file)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(
            hgvs.dataproviders.uta, "connect", lambda *_args, **_kwargs: provider
        )
        # HgvsTools.close() closes its data provider, but this one is shared
        mp.setattr(provider, "close", lambda: None)
        yield provider
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
