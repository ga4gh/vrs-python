import os

import pytest

from ga4gh.vrs.extras.translator import Translator
from biocommons.seqrepo import SeqRepo
from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy, SeqRepoDataProxy
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parents[1]
SEQREPO_ROOT_DIR = os.environ.get("SEQREPO_ROOT_DIR",
                                  f"{APP_ROOT}/tests/data/seqrepo/latest")

@pytest.fixture(scope="session")
def rest_dataproxy():
    return SeqRepoRESTDataProxy(base_url="http://localhost:5000/seqrepo")


@pytest.fixture(scope="session")
def dataproxy():
    sr = SeqRepo(root_dir=SEQREPO_ROOT_DIR)
    return SeqRepoDataProxy(sr)


@pytest.fixture(scope="session")
def rest_tlr(rest_dataproxy):
    return Translator(
        data_proxy=rest_dataproxy,
        default_assembly_name="GRCh38",
        # TODO: Set these to defaults and adjust relevant tests
        identify=False,
        normalize=False,
        translate_sequence_identifiers=True,
    )


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


# See https://github.com/ga4gh/vrs-python/issues/24
# @pytest.fixture(autouse=True)
# def setup_doctest(doctest_namespace, tlr):
#     doctest_namespace["tlr"] = tlr
