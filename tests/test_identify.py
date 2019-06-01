from ga4gh.vr import identify
from ga4gh.vr._internal.computed_identifier import _dictify


class X:
    pass
x = X()
x.id = "an existing id"


def test_use_existing_id():
    "tests that identify() uses an existing id property"
    assert identify(x) == x.id

def test_dictify_None():
    "tests that _dictify works on None"
    assert _dictify(None) == None

def test_dictify_unknown_object():
    "tests that _dictify returns unrecognized objects as-is"
    assert _dictify(x) == x
