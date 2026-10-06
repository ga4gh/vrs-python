import os
from pathlib import Path

import pytest
from biocommons.seqrepo import SeqRepo

from ga4gh.vrs.dataproxy import SeqRepoDataProxy, SeqRepoRESTDataProxy

# Serve hgvs data provider lookups (UTA queries and hgvs sequence fetches) from a
# recorded cache by default, so tests don't need a UTA database or network access. To
# re-record, run `make record-hgvs-cache` with UTA_DB_URL pointing at a UTA instance and
# a seqrepo-rest-service running at SEQREPO_REST_URL.
os.environ.setdefault("VRS_HGVS_CACHE_MODE", "run")
os.environ.setdefault(
    "VRS_HGVS_CACHE_FILE", str(Path(__file__).parent / "data" / "hgvs_cache.pkl")
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
