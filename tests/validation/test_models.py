"""execute models validation tests from the VRS repo

"""

import os

import pytest
import yaml

from ga4gh.core import ga4gh_serialize, ga4gh_digest, ga4gh_identify
from ga4gh.vrs import models

fxs = {
    "ga4gh_serialize": lambda o: ga4gh_serialize(o).decode(),
    "ga4gh_digest": ga4gh_digest,
    "ga4gh_identify": ga4gh_identify,
}

validation_fn = os.path.join(os.path.dirname(__file__), "data", "models.yaml")
validation_tests = yaml.load(open(validation_fn), Loader=yaml.SafeLoader)


def flatten_tests(vts):
    """flatten tests to (class, data, function name, exp out) tuples

    Each tuple is a test in which an object of the specified class is
    created with data, evaluated with the named function, and compared
    with the expected output.

    """
    for cls, tests in validation_tests.items():
        for t in tests:
            for fn, exp in t["out"].items():
                test_name = t.get("name", cls)
                test_name += f"-{fn}"
                yield pytest.param(cls, t["in"], fn, exp, id=test_name)


#tests, ids = zip(*list(flatten_tests(validation_tests)))
#import IPython; IPython.embed()	  ### TODO: Remove IPython.embed()


@pytest.mark.parametrize("cls,data,fn,exp", flatten_tests(validation_tests))
def test_validation(cls, data, fn, exp):
    o = models[cls](**data)
    fx = fxs[fn]
    assert exp == fx(o)
