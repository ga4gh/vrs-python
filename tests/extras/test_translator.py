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
        "sequence": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
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
    "spdi": ["NC_000013.11:20003096:C:", "NC_000013.11:20003096:1:"]
}

deletion_output = {
    "location": {
        "end": 20003097,
        "start": 20003096,
        "sequence": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1687427/?new_evidence=true
insertion_inputs = {
    "hgvs": "NC_000013.11:g.20003010_20003011insG",
    "spdi": ["NC_000013.11:20003010::G", "NC_000013.11:20003010:0:G"]
}

insertion_output = {
    "location": {
        "end": 20003010,
        "start":20003010,
        "sequence": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "G",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1264314/?new_evidence=true
duplication_inputs = {
    "hgvs": "NC_000013.11:g.19993838_19993839dup",
    "spdi": "NC_000013.11:19993837:GT:GTGT"
}

duplication_output = {
    "location": {
        "end": 19993839,
        "start": 19993837,
        "sequence": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "GTGT",
        "type": "LiteralSequenceExpression"
    },
    "type": "Allele"
}


@pytest.mark.vcr
def test_from_beacon(tlr):
    assert tlr._from_beacon(snv_inputs["beacon"]).model_dump(exclude_none=True) == snv_output


@pytest.mark.vcr
def test_from_gnomad(tlr):
    assert tlr._from_gnomad(snv_inputs["gnomad"]).model_dump(exclude_none=True) == snv_output


@pytest.mark.vcr
def test_from_hgvs(tlr):
    assert tlr._from_hgvs(snv_inputs["hgvs"]).model_dump(exclude_none=True) == snv_output
    assert tlr._from_hgvs(deletion_inputs["hgvs"]).model_dump(exclude_none=True) == deletion_output
    assert tlr._from_hgvs(insertion_inputs["hgvs"]).model_dump(exclude_none=True) == insertion_output
    assert tlr._from_hgvs(duplication_inputs["hgvs"]).model_dump(exclude_none=True) == duplication_output


@pytest.mark.vcr
def test_from_spdi(tlr):
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
        "id": "ga4gh:VA.zkZpkDHX0xaOHUPxa5kFt3TUiWbr53xm",
        "location": {
            "id": "ga4gh:SL.2Q2SMLzdWQbJvYbXwwgegrjJbOxiLIJP",
            "end": 32936732,
            "start": 32936731,
            "sequence": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "C",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
    ("NC_000007.14:g.55181320A>T", {
        "id": "ga4gh:VA.RzhTjgnkCmLnaw3IBWnAubZs7eJHhho_",
        "location": {
            "id": "ga4gh:SL.Npx4j5beiNN9GSFTm8Ml6YxrNj_Ghkac",
            "end": 5181320,
            "start": 55181319,
            "sequence": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "T",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
    ("NC_000007.14:g.55181220del", {
        "id": "ga4gh:VA.YdX-N1MsXLxVnfCPYBdkdhzaEJVcPVDA",
        "location": {
            "id": "ga4gh:SL.80HjxOYxb6OtmSUToUehMqU78JJFEaC1",
            "end": 55181220,
            "start": 55181219,
            "sequence": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
    ("NC_000007.14:g.55181230_55181231insGGCT", {
        "id": "ga4gh:VA.osbbI5Me3Paj5mFN3wsD1RpuHhd0dJsv",
        "location": {
            "id": "ga4gh:SL.b30aH2-4AX1oG2LYmMqpBRRyjcpARlgP",
            "end": 55181230,
            "start": 55181230,
            "sequence": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "GGCT",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
    ("NC_000013.11:g.32331093_32331094dup", {
        "id": "ga4gh:VA.0wERl6sqiSnHXe44bfbd4c6zvqAZVjqp",
        "location": {
            "id": "ga4gh:SL.OKUia3LZI-04cFnY1uSO0gEicgTz3suN",
            "end": 32331094,
            "start": 32331082,
            "sequence": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "TTTTTTTTTTTTTT",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
    ("NC_000013.11:g.32316467dup", {
        "id": "ga4gh:VA.iuDk-UaV63YXoYUJHehin3vFfzfYnz0j",
        "location": {
            "id": "ga4gh:SL.yy5lb3KqJGmLbo5kJ7aGf4aW5Ih5-dI3",
            "end": 2316467,
            "start": 32316466,
            "sequence": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "AA",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
    ("NM_001331029.1:n.872A>G", {
        "id": "ga4gh:VA.j-Lnf0txP05DFsLJNhH6Fau4lIBok-DI",
        "location": {
            "id": "ga4gh:SL.VFlsGPHaAJB_LsYzQd8t6vTWSGKWiEeB",
            "end": 872,
            "start": 871,
            "sequence": "ga4gh:SQ.MBIgVnoHFw34aFqNUVGM0zgjC3d-v8dK",
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "G",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }),
    ("NM_181798.1:n.1263G>T", {
        "id": "ga4gh:VA.iLlJQVkxftWXweu4Ad-2NMShYdxUR1DH",
        "location": {
            "id": "ga4gh:SL.rxnYlIIPfVZCxSGhxdvs_nAWlxQd3nen",
            "end": 1263,
            "start": 1262,
            "sequence": "ga4gh:SQ.KN07u-RFqd1dTyOWOG98HnOq87Nq-ZIg",
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
