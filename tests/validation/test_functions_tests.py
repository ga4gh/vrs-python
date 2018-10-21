import os

import pytest
import yaml

import vmc
from vmc.extra.seqrepo import get_vmc_sequence_identifier


validation_fn = os.path.join(os.path.dirname(__file__), "..", "..", "vmc-spec", "tests", "functions.yaml")
validation_tests = yaml.load(open(validation_fn))


@pytest.mark.parametrize("test", validation_tests["digest"])
def test_digest(test):
    assert test["out"]["digest"] == vmc.digest(**test["in"])
