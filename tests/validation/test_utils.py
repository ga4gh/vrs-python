import os

import pytest
import yaml

from ga4gh.vr.extras.seqrepo import get_vmc_sequence_identifier


validation_fn = os.path.join(os.path.dirname(__file__), "..", "data", "schema-tests", "utils.yaml")
validation_tests = yaml.load(open(validation_fn), Loader=yaml.SafeLoader)

@pytest.mark.parametrize("test", validation_tests["get_vmc_sequence_identifier"])
def test_sequence_identifier_translation(test):
    assert test["out"]["identifier"] == get_vmc_sequence_identifier(**test["in"])
