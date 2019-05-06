from ga4gh.vr.extras.translator import Translator

import pytest


@pytest.fixture(scope="session")
def tlr():
    return Translator(default_assembly_name="GRCh38")
