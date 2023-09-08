import pytest

from ga4gh.vrs import models

snv_inputs = {
    "hgvs": "NC_000013.11:g.32936732G>C",
    "beacon": "13 : 32936732 G > C",
    "spdi": "NC_000013.11:32936731:1:C",
    "gnomad": "13-32936732-G-C"
}

snv_output = {
    "location": {
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
        "start":20003010,
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
        "id": "ga4gh:VA.WzOhVdPkrJOnhUL8SftSUstQzkM9oNvJ",
        "location": {
            "id": "ga4gh:SL.4vFSS0UP_RYIDcOihI1qe9D8pfPOCmYO",
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
        "id": "ga4gh:VA.NqnE9Bl89TMfNYExrZbDINwbHlKP9_CT",
        "location": {
            "id": "ga4gh:SL.TQ--CzByOOf6uc89QqNqEzrFLSCkNiv0",
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
        "id": "ga4gh:VA.LkQ5L2nuMTpGHUuX8AQ-ip48R9OoqZyO",
        "location": {
            "id": "ga4gh:SL.UlmvVOCsBly7i3vrmpt9GcHOeyabloZ2",
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
        "id": "ga4gh:VA.wYp6zlK9p7BSwKjSNK4Ntf1t-T7u3lkS",
        "location": {
            "id": "ga4gh:SL.O5Bpm2yq0ca7ytKWMDHM8VxlCxMnWWWR",
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
        "id": "ga4gh:VA.bAl1UhwYzOMjsBRmPTx_n_2vVxLFEB1h",
        "location": {
            "id": "ga4gh:SL.6lpDP0SExMWXuC0-D6TjfPlL4nLgMIqi",
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
        "id": "ga4gh:VA.HOJCOtVRhXGNiRD0Xw1D-QoAU0n5aLi-",
        "location": {
            "id": "ga4gh:SL.ia0FFLAV3CLz-614ow72Dj-xdErUVHg8",
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
        "id": "ga4gh:VA.syC87m7tZ_SE5iIJqfld6wNfp3DdF3aP",
        "location": {
            "id": "ga4gh:SL.fiDihf21EVsPPEoUEG6pTSPbAq7ePxdS",
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
        "id": "ga4gh:VA.rybgbhHiMzGm06azoxsoIkw8SEQ7bdZw",
        "location": {
            "id": "ga4gh:SL.tGsgbZUd5CXglh3KbYV_f6OqxWlIOwqy",
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
    assert expected == allele.model_dump(exclude_none=True)

    to_hgvs = tlr.translate_to(allele, "hgvs")
    assert 1 == len(to_hgvs)
    assert hgvsexpr == to_hgvs[0]


def test_to_hgvs_invalid(tlr):
    # IRI is passed
    iri_vo =  models.Allele(
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
