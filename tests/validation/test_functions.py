from pathlib import Path

import pytest
import yaml

from ga4gh.core import sha512t24u

validation_fn = Path(__file__).parent / "data" / "functions.yaml"
validation_tests = yaml.load(Path(validation_fn).open(), Loader=yaml.SafeLoader)  # noqa: SIM115


@pytest.mark.parametrize("test", validation_tests["sha512t24u"])
def test_sha512t24u(test):
    assert test["out"] == sha512t24u(blob=test["in"]["blob"].encode())
