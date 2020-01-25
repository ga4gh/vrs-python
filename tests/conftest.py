import os

import pytest

from ga4gh.vr.extras.translator import Translator
from ga4gh.vr.extras.dataproxy import SeqRepoRESTDataProxy


@pytest.fixture(scope="session")
def dataproxy():
    return SeqRepoRESTDataProxy(base_url="http://localhost:5000/seqrepo")


@pytest.fixture(scope="session")
def tlr(dataproxy):
    return Translator(
        data_proxy=dataproxy,
        default_assembly_name="GRCh38",
        # TODO: Set these to defaults and adjust relevant tests
        identify=False,
        normalize=False,
        translate_sequence_identifiers=True,
    )


# See https://github.com/ga4gh/vr-python/issues/24
# @pytest.fixture(autouse=True)
# def setup_doctest(doctest_namespace, tlr):
#     doctest_namespace["tlr"] = tlr
