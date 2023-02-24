from copy import deepcopy

import pytest

from ga4gh.core import sha512t24u, ga4gh_digest, ga4gh_serialize, ga4gh_identify, is_pjs_instance
from ga4gh.vrs import models, vrs_deref, vrs_enref

@pytest.fixture(scope="module")
def allele_dict():
    return {
        "type": "Allele",
        "location": {
            "type": "SequenceLocation",
            "interval": {
                "type": "SequenceInterval",
                "start": {
                    "type": "Number",
                    "value": 55181319
                },
                "end": {
                    "type": "Number",
                    "value": 55181320
                }
            },
            "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
        },
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "T"
        }
    }

@pytest.fixture(scope="module")
def allele1():
    a1 = {
        "type": "Allele",
        "location": {
            "type": "SequenceLocation",
            "interval": {
                "type": "SequenceInterval",
                "start": {
                    "type": "Number",
                    "value": 94842865
                },
                "end": {
                    "type": "Number",
                    "value": 94842866
                }
            },
            "sequence_id": "ga4gh:SQ.ss8r_wB0-b9r44TQTMmVTI92884QvBiB"
        },
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "G"
        }
    }
    return models.Allele(**a1)


@pytest.fixture(scope="module")
def allele2():
    a2 = {
        "type": "Allele",
        "location": {
            "type": "SequenceLocation",
            "interval": {
                "type": "SequenceInterval",
                "start": {
                    "type": "Number",
                    "value": 94761899
                },
                "end": {
                    "type": "Number",
                    "value": 94761900
                }
            },
            "sequence_id": "ga4gh:SQ.ss8r_wB0-b9r44TQTMmVTI92884QvBiB"
        },
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "T"
        }
    }
    return models.Allele(**a2)


@pytest.fixture(scope="module")
def incomplete_allele_dict():
    return {
        "type": "Allele",
        "location": {
            "type": "SequenceLocation",
            "interval": {
                "type": "SequenceInterval",
                "start": {
                    "type": "Number",
                    "value": 44908821
                },
                "end": {
                    "type": "Number",
                    "value": 44908822
                }
            },
            "sequence_id": "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
        }
    }


@pytest.fixture(scope="module")
def literal_seq_expr_dict():
    return {
        "type": "LiteralSequenceExpression",
        "sequence": "C"
    }


@pytest.fixture(scope="module")
def repeated_seq_expr_dict():
    return {
        "type": "RepeatedSequenceExpression",
        "count": {
            "type": "IndefiniteRange",
            "value": 6,
            "comparator": ">="
        },
        "seq_expr": {
            "type": "LiteralSequenceExpression",
            "sequence": "A"
        }
    }


def test_vr(allele_dict):
    a = models.Allele(**allele_dict)

    assert a.as_dict() == allele_dict

    assert is_pjs_instance(a.location)

    assert ga4gh_serialize(a.location.interval) == b'{"end":{"type":"Number","value":55181320},"start":{"type":"Number","value":55181319},"type":"SequenceInterval"}'

    assert ga4gh_serialize(
        a.location
    ) == b'{"interval":{"end":{"type":"Number","value":55181320},"start":{"type":"Number","value":55181319},"type":"SequenceInterval"},"sequence_id":"F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceLocation"}'
    assert sha512t24u(ga4gh_serialize(a.location)) == "1tI-5In290Gtp1DOlMxO6z-K3u0OzlXl"
    assert ga4gh_digest(a.location) == "1tI-5In290Gtp1DOlMxO6z-K3u0OzlXl"
    assert ga4gh_identify(a.location) == "ga4gh:VSL.1tI-5In290Gtp1DOlMxO6z-K3u0OzlXl"

    assert ga4gh_serialize(
        a
    ) == b'{"location":"1tI-5In290Gtp1DOlMxO6z-K3u0OzlXl","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}'
    assert ga4gh_digest(a) == "5Z7gWQGUuGAPe4Pw2_kJvnkhS2Q5jRhY"
    assert ga4gh_identify(a) == "ga4gh:VA.5Z7gWQGUuGAPe4Pw2_kJvnkhS2Q5jRhY"

    assert a.as_dict() == allele_dict

    vros = {}
    a2 = vrs_enref(a, vros)
    assert ga4gh_identify(a) == ga4gh_identify(a2)
    assert a2.location == "ga4gh:VSL.1tI-5In290Gtp1DOlMxO6z-K3u0OzlXl"
    assert a2.location in vros
    assert ga4gh_identify(a) in vros

    a3 = vrs_deref(a2, vros)
    assert a == a3
