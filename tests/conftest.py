import os

import pytest
from biocommons.seqrepo import SeqRepo

from ga4gh.vrs.dataproxy import SeqRepoDataProxy, SeqRepoRESTDataProxy


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
