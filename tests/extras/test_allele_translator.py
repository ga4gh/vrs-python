import pytest

from ga4gh.vrs import models
from ga4gh.vrs.extras.translator import AlleleTranslator, ValidationError


@pytest.fixture(scope="module")
def tlr(rest_dataproxy):
    return AlleleTranslator(
        data_proxy=rest_dataproxy,
        default_assembly_name="GRCh38",
        identify=False,
        normalize=False,
    )


snv_inputs = {
    "hgvs": "NC_000019.10:g.44908822C>T",
    "beacon": "19 : 44908822 C > T",
    "spdi": "NC_000019.10:44908821:1:T",
    "gnomad": "19-44908822-C-T"
}

snv_output = {
    "location": {
        "end": 44908822,
        "start": 44908821,
        "sequenceReference": {
            "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "T",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1373966/?new_evidence=true
deletion_inputs = {
    "hgvs": "NC_000013.11:g.20003097del",
    "spdi": ["NC_000013.11:20003096:C:", "NC_000013.11:20003096:1:"],
    "gnomad": "13-20003096-AC-A"
}

deletion_output = {
    "location": {
        "end": 20003097,
        "start": 20003096,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}


gnomad_deletion_output = {
    "location": {
        "end": 20003097,
        "start": 20003095,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "A",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}


deletion_output_normalized = {
    "location": {
        "end": 20003097,
        "start": 20003096,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "length": 0,
        "repeatSubunitLength": 1,
        "sequence": "",
        "type": "ReferenceLengthExpression"
    },
    "type": "Allele"
}


# https://www.ncbi.nlm.nih.gov/clinvar/variation/1687427/?new_evidence=true
insertion_inputs = {
    "hgvs": "NC_000013.11:g.20003010_20003011insG",
    "spdi": ["NC_000013.11:20003010::G", "NC_000013.11:20003010:0:G"],
    "gnomad": "13-20003010-A-AG"
}

insertion_output = {
    "location": {
        "end": 20003010,
        "start": 20003010,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "G",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}


gnomad_insertion_output = {
    "location": {
        "end": 20003010,
        "start": 20003009,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "AG",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}


# https://www.ncbi.nlm.nih.gov/clinvar/variation/1264314/?new_evidence=true
duplication_inputs = {
    "hgvs": "NC_000013.11:g.19993838_19993839dup",
    "spdi": "NC_000013.11:19993837:GT:GTGT",
    "gnomad": "13-19993838-GT-GTGT"
}

duplication_output = {
    "location": {
        "end": 19993839,
        "start": 19993837,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "GTGT",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}


duplication_output_normalized = {
    "location": {
        "end": 19993839,
        "start": 19993837,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "length": 4,
        "repeatSubunitLength": 2,
        "sequence": "GTGT",
        "type": "ReferenceLengthExpression"
    },
    "type": "Allele"
}


@pytest.mark.vcr
def test_from_beacon(tlr):
    tlr.normalize = False
    assert tlr._from_beacon(snv_inputs["beacon"]).model_dump(exclude_none=True) == snv_output


@pytest.mark.vcr
def test_from_gnomad(tlr):
    tlr.normalize = False
    assert tlr._from_gnomad(snv_inputs["gnomad"]).model_dump(exclude_none=True) == snv_output
    assert tlr._from_gnomad(deletion_inputs["gnomad"]).model_dump(exclude_none=True) == gnomad_deletion_output
    assert tlr._from_gnomad(insertion_inputs["gnomad"]).model_dump(exclude_none=True) == gnomad_insertion_output
    assert tlr._from_gnomad(duplication_inputs["gnomad"]).model_dump(exclude_none=True) == duplication_output

    tlr.normalize = True
    assert tlr._from_gnomad(snv_inputs["gnomad"]).model_dump(exclude_none=True) == snv_output
    assert tlr._from_gnomad(deletion_inputs["gnomad"]).model_dump(exclude_none=True) == deletion_output_normalized
    assert tlr._from_gnomad(insertion_inputs["gnomad"]).model_dump(exclude_none=True) == insertion_output
    assert tlr._from_gnomad(duplication_inputs["gnomad"]).model_dump(exclude_none=True) == duplication_output_normalized

    assert tlr._from_gnomad("17-83129587-GTTGWCACATGA-G")

    # Test valid characters
    assert tlr._from_gnomad(
        "7-2-ACGTURYKMSWBDHVN-ACGTURYKMSWBDHVN",
        require_validation=False
    )

    # Invalid input. Ref does not match regex
    assert not tlr._from_gnomad("13-32936732-helloworld-C")

    # Ref != Actual ref
    invalid_var = "13-32936732-G-C"
    error_msg = "Expected reference sequence G on GRCh38:13 at positions (32936731, 32936732) but found C"

    with pytest.raises(ValidationError) as e:
        tlr._from_gnomad(invalid_var)
    assert str(e.value) == error_msg

    with pytest.raises(ValidationError) as e:
        tlr.translate_from(invalid_var, fmt="gnomad")
    assert str(e.value) == error_msg

    # require_validation set to False
    assert tlr._from_gnomad(invalid_var, require_validation=False)


@pytest.mark.vcr
def test_from_hgvs(tlr):
    tlr.normalize = False
    assert tlr._from_hgvs(snv_inputs["hgvs"]).model_dump(exclude_none=True) == snv_output
    assert tlr._from_hgvs(deletion_inputs["hgvs"]).model_dump(exclude_none=True) == deletion_output
    assert tlr._from_hgvs(insertion_inputs["hgvs"]).model_dump(exclude_none=True) == insertion_output
    assert tlr._from_hgvs(duplication_inputs["hgvs"]).model_dump(exclude_none=True) == duplication_output


@pytest.mark.vcr
def test_from_spdi(tlr):
    tlr.normalize = False
    assert tlr._from_spdi(snv_inputs["spdi"]).model_dump(exclude_none=True) == snv_output
    for spdi_del_expr in deletion_inputs["spdi"]:
        assert tlr._from_spdi(spdi_del_expr).model_dump(exclude_none=True) == deletion_output, spdi_del_expr
    for spdi_ins_expr in insertion_inputs["spdi"]:
        assert tlr._from_spdi(spdi_ins_expr).model_dump(exclude_none=True) == insertion_output, spdi_ins_expr
    assert tlr._from_spdi(duplication_inputs["spdi"]).model_dump(exclude_none=True) == duplication_output


@pytest.mark.vcr
def test_to_spdi(tlr):
    tlr.normalize = True
    spdiexpr = snv_inputs["spdi"]
    allele = tlr.translate_from(spdiexpr, "spdi")
    to_spdi = tlr.translate_to(allele, "spdi")
    assert 1 == len(to_spdi)
    assert spdiexpr == to_spdi[0]


hgvs_tests = (
    ("NC_000013.11:g.32936732=", {
        "id": "ga4gh:VA.GuDPEe-WojSx4b4DxupN3si1poaR61qL",
        "location": {
            "id": "ga4gh:SL.jnDT8dINDpAXO35MdLK-8KzLSDL-N3LI",
            "end": 32936732,
            "start": 32936731,
            "sequenceReference": {
                "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "C",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
    ("NC_000007.14:g.55181320A>T", {
        "id": "ga4gh:VA.hOZr7drvRxkUT_srSFVq1NCzvAJdKJlw",
        "location": {
            "id": "ga4gh:SL.wbBBFBTTyBcJPyjkK7z_dCcHFm5pE-2K",
            "end": 55181320,
            "start": 55181319,
            "sequenceReference": {
                "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "T",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
    ("NC_000007.14:g.55181220del", {
        "id": "ga4gh:VA.wlYnlMsWc0ZTPZb-nQv2dXHbFcXa6J9u",
        "location": {
            "id": "ga4gh:SL.hnIOG_kul0Lf3mO1ddTRFb0GbQhtQ19t",
            "end": 55181220,
            "start": 55181219,
            "sequenceReference": {
                "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "state": {
            "length": 0,
            "repeatSubunitLength": 1,
            "sequence": "",
            "type": "ReferenceLengthExpression"
        },
        "type": "Allele"
    }),
    ("NC_000007.14:g.55181230_55181231insGGCT", {
        "id": "ga4gh:VA.lgVw2ZR6UnCLa3sjIMh2D72hQnn6Ksmk",
        "location": {
            "id": "ga4gh:SL.1lu5HsEX4eCQhhoJj_79fjr8H-i98wV3",
            "end": 55181230,
            "start": 55181230,
            "sequenceReference": {
                "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "GGCT",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
    ("NC_000013.11:g.32331093_32331094dup", {
        "id": "ga4gh:VA.x5iNzjjXbb1-wWTBLMBcicYlCMwYoedq",
        "location": {
            "id": "ga4gh:SL.PJ8lHWhAMNRSrxHvkarfDjRWxF-GwaJ_",
            "end": 32331094,
            "start": 32331082,
            "sequenceReference": {
                "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "state": {
            "length": 14,
            "repeatSubunitLength": 2,
            "sequence": "TTTTTTTTTTTTTT",
            "type": "ReferenceLengthExpression"
        },
        "type": "Allele"
    }),
    ("NC_000013.11:g.32316467dup", {
        "id": "ga4gh:VA.ZAyA7Mmd7ERWN6CEd6muxn2mk_gTvEvF",
        "location": {
            "id": "ga4gh:SL.LURTeRdwh5bQf_QqPBoaA--MECYmrY5U",
            "end": 32316467,
            "start": 32316466,
            "sequenceReference": {
                "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "state": {
            "length": 2,
            "repeatSubunitLength": 1,
            "sequence": "AA",
            "type": "ReferenceLengthExpression"
        },
        "type": "Allele"
    }),
    ("NM_001331029.1:n.872A>G", {
        "id": "ga4gh:VA.OKiSoD8f-I0VYwvn9xVGjLxrZ09WWQqK",
        "location": {
            "id": "ga4gh:SL.3flYcxCjrFFW3ex7GqcHLNH8agsKnz49",
            "end": 872,
            "start": 871,
            "sequenceReference": {
                "refgetAccession": "SQ.MBIgVnoHFw34aFqNUVGM0zgjC3d-v8dK",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "G",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
    ("NM_181798.1:n.1263G>T", {
        "id": "ga4gh:VA.wEet5lX0qWJAPuya4MNZsw6ghNYJBSvi",
        "location": {
            "id": "ga4gh:SL.nAdMOa2ccNYPU-DEzzSzN1BaaKvEjYdX",
            "end": 1263,
            "start": 1262,
            "sequenceReference": {
                "refgetAccession": "SQ.KN07u-RFqd1dTyOWOG98HnOq87Nq-ZIg",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "T",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
)


@pytest.mark.parametrize("hgvsexpr,expected", hgvs_tests)
@pytest.mark.vcr
def test_hgvs(tlr, hgvsexpr, expected):
    tlr.normalize = True
    tlr.identify = True
    allele = tlr.translate_from(hgvsexpr, "hgvs")
    assert allele.model_dump(exclude_none=True) == expected

    to_hgvs = tlr.translate_to(allele, "hgvs")
    assert 1 == len(to_hgvs)
    assert hgvsexpr == to_hgvs[0]


def test_to_hgvs_invalid(tlr):
    # IRI is passed
    iri_vo = models.Allele(
        **{
            "location": {
                "end": 1263,
                "start": 1262,
                "sequenceReference": "seqrefs.jsonc#/NM_181798.1",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        }
    )
    with pytest.raises(TypeError) as e:
        tlr.translate_to(iri_vo, "hgvs")
    assert str(e.value) == "`vo.location.sequenceReference` expects a `SequenceReference`"


# TODO: Readd these tests
# @pytest.mark.vcr
# def test_errors(tlr):
#     with pytest.raises(ValueError):
#         tlr._from_beacon("bogus")
#
#     with pytest.raises(ValueError):
#         tlr._from_gnomad("NM_182763.2:c.688+403C>T")
#
#     with pytest.raises(ValueError):
#         tlr._from_hgvs("NM_182763.2:c.688+403C>T")
#
#     with pytest.raises(ValueError):
#         tlr._from_hgvs("NM_182763.2:c.688_690inv")
#
#     with pytest.raises(ValueError):
#         tlr._from_spdi("NM_182763.2:c.688+403C>T")
