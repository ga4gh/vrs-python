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

def test_ordered_false(allele1, allele2):
    """Test that digests are generated correctly when ordered is False"""
    assert ga4gh_serialize(
        allele1.location
    ) == b'{"interval":{"end":{"type":"Number","value":94842866},"start":{"type":"Number","value":94842865},"type":"SequenceInterval"},"sequence_id":"ss8r_wB0-b9r44TQTMmVTI92884QvBiB","type":"SequenceLocation"}'
    assert sha512t24u(ga4gh_serialize(allele1.location)) == "yNuwN3VJz77KdF6mY0GN_S5bnw7ByV2L"
    assert ga4gh_digest(allele1.location) == "yNuwN3VJz77KdF6mY0GN_S5bnw7ByV2L"
    assert ga4gh_identify(allele1.location) == "ga4gh:VSL.yNuwN3VJz77KdF6mY0GN_S5bnw7ByV2L"

    assert ga4gh_serialize(
        allele1
    ) == b'{"location":"yNuwN3VJz77KdF6mY0GN_S5bnw7ByV2L","state":{"sequence":"G","type":"LiteralSequenceExpression"},"type":"Allele"}'
    assert ga4gh_digest(allele1) == "geQCxa1Enel8UBUAQQ2-rbphDjIR-cq0"
    assert ga4gh_identify(allele1) == "ga4gh:VA.geQCxa1Enel8UBUAQQ2-rbphDjIR-cq0"

    assert ga4gh_serialize(
        allele2.location
    ) == b'{"interval":{"end":{"type":"Number","value":94761900},"start":{"type":"Number","value":94761899},"type":"SequenceInterval"},"sequence_id":"ss8r_wB0-b9r44TQTMmVTI92884QvBiB","type":"SequenceLocation"}'
    assert sha512t24u(ga4gh_serialize(allele2.location)) == "2b33Z0JsziOXLGOtKOu-jPwz0VL9Gfp-"
    assert ga4gh_digest(allele2.location) == "2b33Z0JsziOXLGOtKOu-jPwz0VL9Gfp-"
    assert ga4gh_identify(allele2.location) == "ga4gh:VSL.2b33Z0JsziOXLGOtKOu-jPwz0VL9Gfp-"

    assert ga4gh_serialize(
        allele2
    ) == b'{"location":"2b33Z0JsziOXLGOtKOu-jPwz0VL9Gfp-","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}'
    assert ga4gh_digest(allele2) == "jWqv036CdZJs4YjwEYptDIBcoT7Uxv5I"
    assert ga4gh_identify(allele2) == "ga4gh:VA.jWqv036CdZJs4YjwEYptDIBcoT7Uxv5I"

    expected_haplotype_serialization = b'{"members":["geQCxa1Enel8UBUAQQ2-rbphDjIR-cq0","jWqv036CdZJs4YjwEYptDIBcoT7Uxv5I"],"type":"Haplotype"}'
    expected_haplotype_digest = "Ow_uE0YaVWHIno4pQfdmYpWmlGPNtXQr"
    expected_haplotype_curie = f"ga4gh:VH.{expected_haplotype_digest}"

    h1 = models.Haplotype(
        type="Haplotype",
        members=[allele1, allele2]
    )
    assert ga4gh_serialize(h1) == expected_haplotype_serialization
    assert ga4gh_digest(h1) == expected_haplotype_digest
    assert ga4gh_identify(h1) == expected_haplotype_curie

    h2 = models.Haplotype(
        type="Haplotype",
        members=[allele2, allele1]
    )
    assert ga4gh_serialize(h2) == expected_haplotype_serialization
    assert ga4gh_digest(h2) == expected_haplotype_digest
    assert ga4gh_identify(h2) == expected_haplotype_curie

    assert ga4gh_identify(h1) == ga4gh_identify(h2) == expected_haplotype_curie

def test_ordered_true(incomplete_allele_dict, literal_seq_expr_dict, repeated_seq_expr_dict):
    """Test that digests are generated correctly when ordered is True"""
    a1 = deepcopy(incomplete_allele_dict)
    a1["state"] = {
        "type": "ComposedSequenceExpression",
        "components": [literal_seq_expr_dict, repeated_seq_expr_dict]
    }
    a1 = models.Allele(**a1)
    assert ga4gh_serialize(
        a1.location
    ) == b'{"interval":{"end":{"type":"Number","value":44908822},"start":{"type":"Number","value":44908821},"type":"SequenceInterval"},"sequence_id":"IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl","type":"SequenceLocation"}'
    assert sha512t24u(ga4gh_serialize(a1.location)) == "QrRSuBj-VScAGV_gEdxNgsnh41jYH1Kg"
    assert ga4gh_digest(a1.location) == "QrRSuBj-VScAGV_gEdxNgsnh41jYH1Kg"
    assert ga4gh_identify(a1.location) == "ga4gh:VSL.QrRSuBj-VScAGV_gEdxNgsnh41jYH1Kg"

    assert ga4gh_serialize(
        a1
    ) == b'{"location":"QrRSuBj-VScAGV_gEdxNgsnh41jYH1Kg","state":{"components":[{"sequence":"C","type":"LiteralSequenceExpression"},{"count":{"comparator":">=","type":"IndefiniteRange","value":6},"seq_expr":{"sequence":"A","type":"LiteralSequenceExpression"},"type":"RepeatedSequenceExpression"}],"type":"ComposedSequenceExpression"},"type":"Allele"}'
    assert ga4gh_digest(a1) == "HQeEhxKyJNirJAUlgyq-psa7UlE9l9X3"
    assert ga4gh_identify(a1) == "ga4gh:VA.HQeEhxKyJNirJAUlgyq-psa7UlE9l9X3"


    a2 = deepcopy(incomplete_allele_dict)
    a2["state"] = {
        "type": "ComposedSequenceExpression",
        "components": [repeated_seq_expr_dict, literal_seq_expr_dict]
    }
    a2 = models.Allele(**a2)
    assert ga4gh_identify(a2.location) == ga4gh_identify(a1.location)

    assert ga4gh_serialize(
        a2
    ) == b'{"location":"QrRSuBj-VScAGV_gEdxNgsnh41jYH1Kg","state":{"components":[{"count":{"comparator":">=","type":"IndefiniteRange","value":6},"seq_expr":{"sequence":"A","type":"LiteralSequenceExpression"},"type":"RepeatedSequenceExpression"},{"sequence":"C","type":"LiteralSequenceExpression"}],"type":"ComposedSequenceExpression"},"type":"Allele"}'
    assert ga4gh_digest(a2) == "_rPRyTJwcgBmSrdtgrEcOmRY9qZLkc5N"
    assert ga4gh_identify(a2) == "ga4gh:VA._rPRyTJwcgBmSrdtgrEcOmRY9qZLkc5N"

    assert ga4gh_identify(a1) != ga4gh_identify(a2)
