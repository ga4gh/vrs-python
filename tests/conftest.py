import os

import pytest

from ga4gh.vrs.extras.translator import Translator
from biocommons.seqrepo import SeqRepo
from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy, SeqRepoDataProxy


@pytest.fixture(scope="session")
def dataproxy():
    sr = SeqRepo(root_dir=os.environ.get("SEQREPO_ROOT_DIR", "/usr/local/share/seqrepo/latest"))
    return SeqRepoDataProxy(sr)


@pytest.fixture(scope="session")
def rest_dataproxy():
    return SeqRepoRESTDataProxy(base_url="http://localhost:5000/seqrepo")


@pytest.fixture(scope="session")
def tlr(rest_dataproxy):
    return Translator(
        data_proxy=rest_dataproxy,
        default_assembly_name="GRCh38",
    # TODO: Set these to defaults and adjust relevant tests
        identify=False,
        normalize=False,
        translate_sequence_identifiers=True,
    )


# See https://github.com/ga4gh/vrs-python/issues/24
# @pytest.fixture(autouse=True)
# def setup_doctest(doctest_namespace, tlr):
#     doctest_namespace["tlr"] = tlr
