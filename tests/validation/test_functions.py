import os

import pytest
import vcr
import yaml

from ga4gh.core import sha512t24u

validation_fn = os.path.join(os.path.dirname(__file__), "data", "functions.yaml")
validation_tests = yaml.load(open(validation_fn), Loader=yaml.SafeLoader)


@pytest.mark.parametrize("test", validation_tests["sha512t24u"])
def test_sha512t24u(test):
    assert test["out"] == sha512t24u(blob=test["in"]["blob"].encode())
