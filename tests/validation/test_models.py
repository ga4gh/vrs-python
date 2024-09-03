"""execute models validation tests from the VRS repo

"""

import os

from pydantic import ValidationError
import pytest
import yaml

from ga4gh.core import ga4gh_serialize, ga4gh_digest, ga4gh_identify, PrevVrsVersion, entity_models, CommonDomainType, domain_models
from ga4gh.vrs import models, VrsType

def ga4gh_1_3_identify(*args, **kwargs):
    kwargs['as_version'] = PrevVrsVersion.V1_3
    return ga4gh_identify(*args, **kwargs)

def ga4gh_1_3_digest(*args, **kwargs):
    kwargs['as_version'] = PrevVrsVersion.V1_3
    return ga4gh_digest(*args, **kwargs)

def ga4gh_1_3_serialize(*args, **kwargs):
    kwargs['as_version'] = PrevVrsVersion.V1_3
    return ga4gh_serialize(*args, **kwargs)

fxs = {
    "ga4gh_serialize": ga4gh_serialize,
    "ga4gh_digest": ga4gh_digest,
    "ga4gh_identify": ga4gh_identify,
    "ga4gh_1_3_digest": ga4gh_1_3_digest,
    "ga4gh_1_3_identify": ga4gh_1_3_identify,
    "ga4gh_1_3_serialize": ga4gh_1_3_serialize
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
    o = getattr(models, cls)(**data)
    fx = fxs[fn]
    if fn == "ga4gh_serialize":
        exp = exp.encode("utf-8")
    assert fx(o) == exp


def test_prev_vrs_version():
    """Ensure that support to previous VRS digest/identifiers works correctly"""
    loc = models.SequenceLocation(start=44908821, end=44908822, sequenceReference=models.SequenceReference(refgetAccession="SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"))

    # string representation should work as well
    ga4gh_identify(loc, as_version="1.3")

    invalid_vrs_version = "0.0"
    invalid_vrs_version_msg = f"Expected `PrevVrsVersion`, but got {invalid_vrs_version}"

    loc_no_seq_ref = models.SequenceLocation(start=44908821, end=44908822)
    loc_iri = models.SequenceLocation(start=44908821, end=44908822, sequenceReference=entity_models.IRI("sequenceReferences.json#example1"))
    allele_rle_no_seq = models.Allele(location=loc, state=models.ReferenceLengthExpression(length=11, repeatSubunitLength=3))
    allele_le = models.Allele(location=loc, state=models.LengthExpression(length=2))
    loc_seq_ref_msg = "Must provide `sequenceReference` and it must be a valid `SequenceReference`"
    for ga4gh_func in [ga4gh_identify, ga4gh_digest, ga4gh_serialize]:
        with pytest.raises(ValueError, match=invalid_vrs_version_msg):
            ga4gh_func(loc, as_version=invalid_vrs_version_msg)

        with pytest.raises(ValueError, match=loc_seq_ref_msg):
            ga4gh_func(loc_no_seq_ref, as_version=PrevVrsVersion.V1_3)

        with pytest.raises(ValueError, match=loc_seq_ref_msg):
            ga4gh_func(loc_iri, as_version=PrevVrsVersion.V1_3)

        with pytest.raises(ValueError, match="State sequence attribute must be defined."):
            ga4gh_func(allele_rle_no_seq, as_version=PrevVrsVersion.V1_3)

        allele_rlse_seq = allele_rle_no_seq.model_copy(deep=True)
        allele_rlse_seq.state.sequence = "C"
        assert ga4gh_func(allele_rlse_seq, as_version=PrevVrsVersion.V1_3)

        with pytest.raises(ValueError, match="Only `LiteralSequenceExpression` and `ReferenceLengthExpression` are supported for previous versions of VRS"):
            ga4gh_func(allele_le, as_version=PrevVrsVersion.V1_3)


def test_valid_types():
    """Ensure that type enums values correct. Values should correspond to class"""
    for gks_models, gks_enum in [(models, VrsType), (domain_models, CommonDomainType)]:
        for enum_val in gks_enum.__members__.values():
            enum_val = enum_val.value
            if hasattr(gks_models, enum_val):
                gks_class = getattr(gks_models, enum_val)
                try:
                    assert gks_class(type=enum_val)
                except ValidationError as e:
                    found_type_mismatch = False
                    for error in e.errors():
                        if error["loc"] == ("type",):
                            found_type_mismatch = True
                    assert not found_type_mismatch, f"Found mismatch in type literal: {enum_val} vs {error['ctx']['expected']}"
            else:
                assert False, f"{str(gks_models)} class not found: {enum_val}"
