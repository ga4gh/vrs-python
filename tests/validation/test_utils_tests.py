import os

import pytest
import yaml

import vmc

validation_fn = os.path.join(os.path.dirname(__file__), "..", "..", "vmc-spec", "tests", "utils.yaml")
validation_tests = yaml.load(open(validation_fn), Loader=yaml.SafeLoader)


@pytest.mark.parametrize("test", validation_tests["get_vmc_sequence_identifier"])
def test_sequence_identifier_translation(test):
    assert test["out"]["identifier"] == vmc.get_vmc_sequence_identifier(**test["in"])
