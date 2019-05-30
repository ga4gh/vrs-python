import os

import pytest
import vcr

from ga4gh.vr.extras.translator import Translator
from ga4gh.vr.extras.dataproxy import SeqRepoRESTDataProxy


test_dir = os.path.dirname(__file__)
test_data_dir = os.path.join(test_dir, "data", "cassettes")

vcr.default_vcr = vcr.VCR(
    cassette_library_dir=test_data_dir,
    record_mode=os.environ.get("VCR_RECORD_MODE", "once"),
    )
vcr.use_cassette = vcr.default_vcr.use_cassette


@vcr.use_cassette
@pytest.fixture(scope="session")
def dataproxy():
    return SeqRepoRESTDataProxy(base_url="http://localhost:5000/seqrepo")


@pytest.fixture(scope="session")
def tlr(dataproxy):
    return Translator(data_proxy=dataproxy,
                      default_assembly_name="GRCh38")

