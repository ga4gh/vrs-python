import os

import pytest

from biocommons.seqrepo import SeqRepo
from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy, SeqRepoDataProxy


@pytest.fixture(scope="session")
def dataproxy():
    sr = SeqRepo(root_dir=os.environ.get("SEQREPO_ROOT_DIR", "/usr/local/share/seqrepo/latest"))
    return SeqRepoDataProxy(sr)


@pytest.fixture(scope="session")
def rest_dataproxy():
    return SeqRepoRESTDataProxy(
        base_url=os.environ.get(
            "SEQREPO_REST_URL",
            "http://localhost:5000/seqrepo"))
