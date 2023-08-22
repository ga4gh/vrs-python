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
        "id": "ga4gh:VA.Xj-QkZI6fXYy6vW0XejUBHl1cfdCF4m3",
        "location": {
            "id": "ga4gh:SL.hRWzC4ss2ici8DtHbdvj-1ApbDFRK2Yd",
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
        "id": "ga4gh:VA.Qifv6-W3fIeg0_yn53tu2p-MC4wnLPXr",
        "location": {
            "id": "ga4gh:SL.7qnnS0w7J8d2-iEwuP6qxZYUBVp07HDO",
            "end": 55181320,
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
        "id": "ga4gh:VA.WLf4uNKhsBMp2A6Xm43ZzkPzmXYFPiRd",
        "location": {
            "id": "ga4gh:SL.EPAFk0DH8VQ8ItAOTyUrMXHQAAOEO4SV",
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
        "id": "ga4gh:VA.KMg70UdX0kanxa2fj7iKvL3YNtVOGk6_",
        "location": {
            "id": "ga4gh:SL.zqIGqMt8l8q2LQrUAzu4IR9Xe_72Z8lV",
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
        "id": "ga4gh:VA.1QW7zxUdGhQbr2AQ0QKO5d329_SQCWz-",
        "location": {
            "id": "ga4gh:SL.1pzga0d1li6jaFVeE_H9BndUHcpxxO4n",
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
        "id": "ga4gh:VA.eQzjNd7OOgIaB7V5U8HL9njU6hOB_UeH",
        "location": {
            "id": "ga4gh:SL.yERDq_Lr1yCMiuM75ALuCuZBpsI-y4wN",
            "end": 32316467,
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
        "id": "ga4gh:VA.8sRg3IThjwkZmxjYMxiyMrbiA2H6L51k",
        "location": {
            "id": "ga4gh:SL.YunNaSFmYIQAVrrYHaoC9lMT4FAcUovg",
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
        "id": "ga4gh:VA.fs1txJnlahHtGCSgiVRgHW9bLg0DNaAI",
        "location": {
            "id": "ga4gh:SL.pzDOS6YY0BLGeTE9HlL1d66uxBc5Hd6l",
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
