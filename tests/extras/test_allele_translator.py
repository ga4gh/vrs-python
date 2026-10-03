from typing import ClassVar
from unittest.mock import MagicMock, patch

import pytest

from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import DataProxyValidationError, _DataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.utils.hgvs_tools import HgvsTools


@pytest.fixture(scope="module")
def tlr(rest_dataproxy):
    return AlleleTranslator(
        data_proxy=rest_dataproxy,
        default_assembly_name="GRCh38",
        identify=False,
    )


# https://www.ncbi.nlm.nih.gov/clinvar/variation/17848/?new_evidence=true
snv_inputs = {
    "hgvs": "NC_000019.10:g.44908822C>T",
    "beacon": "19 : 44908822 C > T",
    "spdi": "NC_000019.10:44908821:1:T",
    "gnomad": "19-44908822-C-T",
}

snv_output = {
    "location": {
        "end": 44908822,
        "start": 44908821,
        "sequenceReference": {
            "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
            "type": "SequenceReference",
        },
        "type": "SequenceLocation",
    },
    "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
    "type": "Allele",
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/693259/?new_evidence=true
mito_inputs = {
    "hgvs": "NC_012920.1:m.10083A>G",
    "beacon": "MT : 10083 A > G",
    "spdi": "NC_012920.1:10082:A:G",
    "gnomad": "MT-10083-A-G",
}

mito_output = {
    "location": {
        "start": 10082,
        "end": 10083,
        "sequenceReference": {
            "refgetAccession": "SQ.k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct",
            "type": "SequenceReference",
        },
        "type": "SequenceLocation",
    },
    "state": {"sequence": "G", "type": "LiteralSequenceExpression"},
    "type": "Allele",
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1373966/?new_evidence=true
deletion_inputs = {
    "hgvs": "NC_000013.11:g.20003097del",
    "spdi": ["NC_000013.11:20003096:C:", "NC_000013.11:20003096:1:"],
    "gnomad": "13-20003096-AC-A",
}

deletion_output = {
    "location": {
        "end": 20003097,
        "start": 20003096,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference",
        },
        "type": "SequenceLocation",
    },
    "state": {"sequence": "", "type": "LiteralSequenceExpression"},
    "type": "Allele",
}

gnomad_deletion_output = {
    "location": {
        "end": 20003097,
        "start": 20003095,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference",
        },
        "type": "SequenceLocation",
    },
    "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
    "type": "Allele",
}

deletion_output_normalized = {
    "location": {
        "end": 20003097,
        "start": 20003096,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference",
        },
        "type": "SequenceLocation",
    },
    "state": {
        "length": 0,
        "repeatSubunitLength": 1,
        "sequence": "",
        "type": "ReferenceLengthExpression",
    },
    "type": "Allele",
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1687427/?new_evidence=true
insertion_inputs = {
    "hgvs": "NC_000013.11:g.20003010_20003011insG",
    "spdi": ["NC_000013.11:20003010::G", "NC_000013.11:20003010:0:G"],
    "gnomad": "13-20003010-A-AG",
}

insertion_output = {
    "location": {
        "end": 20003010,
        "start": 20003010,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference",
        },
        "type": "SequenceLocation",
    },
    "state": {"sequence": "G", "type": "LiteralSequenceExpression"},
    "type": "Allele",
}

gnomad_insertion_output = {
    "location": {
        "end": 20003010,
        "start": 20003009,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference",
        },
        "type": "SequenceLocation",
    },
    "state": {"sequence": "AG", "type": "LiteralSequenceExpression"},
    "type": "Allele",
}

# https://www.ncbi.nlm.nih.gov/clinvar/variation/1264314/?new_evidence=true
duplication_inputs = {
    "hgvs": "NC_000013.11:g.19993838_19993839dup",
    "spdi": "NC_000013.11:19993837:GT:GTGT",
    "gnomad": "13-19993838-GT-GTGT",
}

duplication_output = {
    "location": {
        "end": 19993839,
        "start": 19993837,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference",
        },
        "type": "SequenceLocation",
    },
    "state": {"sequence": "GTGT", "type": "LiteralSequenceExpression"},
    "type": "Allele",
}

duplication_output_normalized = {
    "location": {
        "end": 19993839,
        "start": 19993837,
        "sequenceReference": {
            "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
            "type": "SequenceReference",
        },
        "type": "SequenceLocation",
    },
    "state": {
        "length": 4,
        "repeatSubunitLength": 2,
        "sequence": "GTGT",
        "type": "ReferenceLengthExpression",
    },
    "type": "Allele",
}

rle_inputs = [
    {
        # small 1+ repeat deletion
        "gnomad": "1-145916840-CTCCT-CT",
        "spdi": "NC_000001.11:145916839:CTCCT:CT",
        "expected_allele": {
            "type": "Allele",
            "location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                },
                "start": 145916839,
                "end": 145916844,
            },
            "state": {
                "type": "ReferenceLengthExpression",
                "length": 2,
                "sequence": "CT",
                "repeatSubunitLength": 3,
            },
        },
    },
    {
        # big 2+ multiple-repeat deletion
        "gnomad": "21-5033800-TGGGTCCAGGCACCGGCGCCCAGCCCCCGTGGGGTGTCCAGGGC-T",
        "spdi": "NC_000021.9:5033800:GGGTCCAGGCACCGGCGCCCAGCCCCCGTGGGGTGTCCAGGGCGGGTCCAGGCACCGGCGCCCAGCCCCCGTGGGGTGTCCAGGGCGGGTCCAGGCACCGGCGCCCAGCCCCC:GGGTCCAGGCACCGGCGCCCAGCCCCCGTGGGGTGTCCAGGGCGGGTCCAGGCACCGGCGCCCAGCCCCC",
        "expected_allele": {
            "location": {
                "start": 5033800,
                "end": 5033913,
                "sequenceReference": {
                    "refgetAccession": "SQ.5ZUqxCmDDgN4xTRbaSjN8LwgZironmB8",
                    "type": "SequenceReference",
                },
                "type": "SequenceLocation",
            },
            "state": {
                "length": 70,
                "repeatSubunitLength": 43,
                "type": "ReferenceLengthExpression",
                "sequence": "GGGTCCAGGCACCGGCGCCCAGCCCCCGTGGGGTGTCCAGGGCGGGTCCAGGCACCGGCGCCCAGCCCCC",
            },
            "type": "Allele",
        },
    },
    {
        # insertion
        "gnomad": "1-2228956-GTGCCCG-GTGCCCGTGCCCG",
        "spdi": "NC_000001.11:2228955:GTGCCCG:GTGCCCGTGCCCG",
        "expected_allele": {
            "type": "Allele",
            "location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                },
                "start": 2228955,
                "end": 2228962,
            },
            "state": {
                "type": "ReferenceLengthExpression",
                "length": 13,
                "sequence": "GTGCCCGTGCCCG",
                "repeatSubunitLength": 6,
            },
        },
    },
    {
        # single base expansion (9*A -> 11*A)
        "spdi": "NC_000001.11:236900409:AAAAAAAAA:AAAAAAAAAAA",
        "gnomad": "1-236900410-A-AAA",
        "expected_allele": {
            "location": {
                "end": 236900418,
                "sequenceReference": {
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    "type": "SequenceReference",
                },
                "start": 236900409,
                "type": "SequenceLocation",
            },
            "state": {
                "length": 11,
                "repeatSubunitLength": 2,
                "sequence": "AAAAAAAAAAA",
                "type": "ReferenceLengthExpression",
            },
            "type": "Allele",
        },
    },
    {
        # single base expansion (9*A -> 11*A)
        "spdi": "NC_000001.11:236900409:AAAAAAAAA:AAAAAAAAAAAAAAAAAAAA",
        "gnomad": "1-236900410-A-AAAAAAAAAAAA",
        "expected_allele": {
            "location": {
                "end": 236900418,
                "sequenceReference": {
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    "type": "SequenceReference",
                },
                "start": 236900409,
                "type": "SequenceLocation",
            },
            "state": {
                "length": 20,
                "repeatSubunitLength": 1,
                "sequence": "AAAAAAAAAAAAAAAAAAAA",
                "type": "ReferenceLengthExpression",
            },
            "type": "Allele",
        },
    },
]


@pytest.mark.vcr
def test_rle_round_trip_gnomad_spdi(tlr):
    """Test translating an RLE allele from gnomAD and SPDI, and back to SPDI"""
    for inp in rle_inputs:
        gnomad = inp["gnomad"]
        spdi = inp["spdi"]
        expected_dict = inp["expected_allele"]
        allele_gnomad = tlr.translate_from(gnomad, fmt="gnomad", rle_seq_limit=None)
        assert allele_gnomad.model_dump(exclude_none=True) == expected_dict
        allele_spdi = tlr.translate_from(spdi, fmt="spdi", rle_seq_limit=None)
        assert allele_spdi.model_dump(exclude_none=True) == expected_dict

        allele_gnomad_to_spdi = tlr.translate_to(
            allele_gnomad, fmt="spdi", ref_seq_limit=None
        )
        assert len(allele_gnomad_to_spdi) == 1
        assert allele_gnomad_to_spdi[0] == spdi


def test_from_invalid(tlr):
    with pytest.raises(
        ValueError, match="Unable to parse data as beacon, gnomad, hgvs, spdi, vrs"
    ):
        tlr.translate_from("BRAF amplication")

    with pytest.raises(
        ValueError, match="Unable to parse data as beacon, gnomad, hgvs, spdi, vrs"
    ):
        tlr.translate_from("BRAF amplication", assembly_name="GRCh37")


@pytest.mark.vcr
def test_from_beacon(tlr):
    do_normalize = False
    with pytest.deprecated_call():
        assert (
            tlr._from_beacon(
                snv_inputs["beacon"], do_normalize=do_normalize
            ).model_dump(exclude_none=True)
            == snv_output
        )

    with pytest.deprecated_call():
        assert (
            tlr._from_beacon(
                mito_inputs["beacon"], do_normalize=do_normalize
            ).model_dump(exclude_none=True)
            == mito_output
        )


@pytest.mark.vcr
def test_from_gnomad(tlr):
    do_normalize = False
    assert (
        tlr._from_gnomad(snv_inputs["gnomad"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == snv_output
    )
    assert (
        tlr._from_gnomad(mito_inputs["gnomad"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == mito_output
    )
    assert (
        tlr._from_gnomad(
            deletion_inputs["gnomad"], do_normalize=do_normalize
        ).model_dump(exclude_none=True)
        == gnomad_deletion_output
    )
    assert (
        tlr._from_gnomad(
            insertion_inputs["gnomad"], do_normalize=do_normalize
        ).model_dump(exclude_none=True)
        == gnomad_insertion_output
    )
    assert (
        tlr._from_gnomad(
            duplication_inputs["gnomad"], do_normalize=do_normalize
        ).model_dump(exclude_none=True)
        == duplication_output
    )

    # do_normalize defaults to true
    assert (
        tlr._from_gnomad(snv_inputs["gnomad"]).model_dump(exclude_none=True)
        == snv_output
    )
    assert (
        tlr._from_gnomad(mito_inputs["gnomad"]).model_dump(exclude_none=True)
        == mito_output
    )
    assert (
        tlr._from_gnomad(deletion_inputs["gnomad"]).model_dump(exclude_none=True)
        == deletion_output_normalized
    )
    assert (
        tlr._from_gnomad(insertion_inputs["gnomad"]).model_dump(exclude_none=True)
        == insertion_output
    )
    assert (
        tlr._from_gnomad(duplication_inputs["gnomad"]).model_dump(exclude_none=True)
        == duplication_output_normalized
    )

    assert tlr._from_gnomad("17-83129587-GTTGWCACATGA-G")

    # Test valid characters
    assert tlr._from_gnomad(
        "7-2-ACGTURYKMSWBDHVN-ACGTURYKMSWBDHVN", require_validation=False
    )

    # Invalid input. Ref does not match regex
    assert not tlr._from_gnomad("13-32936732-helloworld-C")

    # Ref != Actual ref
    invalid_var = "13-32936732-G-C"
    error_msg = "Reference mismatch at GRCh38:13 position 32936731-32936732 (input gave 'G' but correct ref is 'C')"

    with pytest.raises(DataProxyValidationError) as e:
        tlr._from_gnomad(invalid_var)
    assert str(e.value) == error_msg

    with pytest.raises(DataProxyValidationError) as e:
        tlr.translate_from(invalid_var, fmt="gnomad")
    assert str(e.value) == error_msg

    # require_validation set to False
    assert tlr._from_gnomad(invalid_var, require_validation=False)


@pytest.mark.vcr
def test_from_hgvs(tlr):
    do_normalize = False
    assert (
        tlr._from_hgvs(snv_inputs["hgvs"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == snv_output
    )
    assert (
        tlr._from_hgvs(mito_inputs["hgvs"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == mito_output
    )
    assert (
        tlr._from_hgvs(deletion_inputs["hgvs"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == deletion_output
    )
    assert (
        tlr._from_hgvs(insertion_inputs["hgvs"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == insertion_output
    )
    assert (
        tlr._from_hgvs(
            duplication_inputs["hgvs"], do_normalize=do_normalize
        ).model_dump(exclude_none=True)
        == duplication_output
    )


@pytest.mark.vcr
def test_from_spdi(tlr):
    do_normalize = False
    assert (
        tlr._from_spdi(snv_inputs["spdi"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == snv_output
    )
    assert (
        tlr._from_spdi(mito_inputs["spdi"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == mito_output
    )
    for spdi_del_expr in deletion_inputs["spdi"]:
        assert (
            tlr._from_spdi(spdi_del_expr, do_normalize=do_normalize).model_dump(
                exclude_none=True
            )
            == deletion_output
        ), spdi_del_expr
    for spdi_ins_expr in insertion_inputs["spdi"]:
        assert (
            tlr._from_spdi(spdi_ins_expr, do_normalize=do_normalize).model_dump(
                exclude_none=True
            )
            == insertion_output
        ), spdi_ins_expr
    assert (
        tlr._from_spdi(
            duplication_inputs["spdi"], do_normalize=do_normalize
        ).model_dump(exclude_none=True)
        == duplication_output
    )


@pytest.mark.vcr
def test_to_spdi(tlr):
    # do_normalize defaults to true
    spdiexpr = snv_inputs["spdi"]
    allele = tlr.translate_from(spdiexpr, "spdi")
    to_spdi = tlr.translate_to(allele, "spdi")
    assert len(to_spdi) == 1
    assert spdiexpr == to_spdi[0]


@pytest.mark.vcr
def test_to_spdi_with_ref(tlr):
    spdi_expr_no_ref = "NC_000019.10:44908821:1:T"
    spdi_expr_with_ref = "NC_000019.10:44908821:C:T"

    allele_no_ref = tlr.translate_from(spdi_expr_no_ref, "spdi")
    allele_with_ref = tlr.translate_from(spdi_expr_with_ref, "spdi")
    assert allele_no_ref == allele_with_ref

    to_spdi_no_ref = tlr.translate_to(allele_no_ref, "spdi", ref_seq_limit=0)
    assert len(to_spdi_no_ref) == 1
    assert spdi_expr_no_ref == to_spdi_no_ref[0]

    to_spdi_with_ref = tlr.translate_to(allele_with_ref, "spdi", ref_seq_limit=None)
    assert len(to_spdi_with_ref) == 1
    assert spdi_expr_with_ref == to_spdi_with_ref[0]


hgvs_tests = (
    (
        "NC_000013.11:g.32936732=",
        {
            "location": {
                "end": 32936732,
                "sequenceReference": {
                    "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                    "type": "SequenceReference",
                },
                "start": 32936731,
                "type": "SequenceLocation",
            },
            "state": {
                "length": 1,
                "repeatSubunitLength": 1,
                "sequence": "C",
                "type": "ReferenceLengthExpression",
            },
            "type": "Allele",
        },
    ),
    (
        "NC_000007.14:g.55181320A>T",
        {
            "location": {
                "end": 55181320,
                "sequenceReference": {
                    "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "type": "SequenceReference",
                },
                "start": 55181319,
                "type": "SequenceLocation",
            },
            "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
    ),
    (
        "NC_000007.14:g.55181220del",
        {
            "location": {
                "end": 55181220,
                "sequenceReference": {
                    "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "type": "SequenceReference",
                },
                "start": 55181219,
                "type": "SequenceLocation",
            },
            "state": {
                "length": 0,
                "repeatSubunitLength": 1,
                "sequence": "",
                "type": "ReferenceLengthExpression",
            },
            "type": "Allele",
        },
    ),
    (
        "NC_000007.14:g.55181230_55181231insGGCT",
        {
            "location": {
                "end": 55181230,
                "sequenceReference": {
                    "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "type": "SequenceReference",
                },
                "start": 55181230,
                "type": "SequenceLocation",
            },
            "state": {"sequence": "GGCT", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
    ),
    (
        "NC_000013.11:g.32331093_32331094dup",
        {
            "location": {
                "end": 32331094,
                "sequenceReference": {
                    "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                    "type": "SequenceReference",
                },
                "start": 32331082,
                "type": "SequenceLocation",
            },
            "state": {
                "length": 14,
                "repeatSubunitLength": 2,
                "sequence": "TTTTTTTTTTTTTT",
                "type": "ReferenceLengthExpression",
            },
            "type": "Allele",
        },
    ),
    (
        "NC_000013.11:g.32316467dup",
        {
            "location": {
                "end": 32316467,
                "sequenceReference": {
                    "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                    "type": "SequenceReference",
                },
                "start": 32316466,
                "type": "SequenceLocation",
            },
            "state": {
                "length": 2,
                "repeatSubunitLength": 1,
                "sequence": "AA",
                "type": "ReferenceLengthExpression",
            },
            "type": "Allele",
        },
    ),
    (
        "NM_001331029.1:c.722A>G",
        {
            "location": {
                "end": 872,
                "sequenceReference": {
                    "refgetAccession": "SQ.MBIgVnoHFw34aFqNUVGM0zgjC3d-v8dK",
                    "type": "SequenceReference",
                },
                "start": 871,
                "type": "SequenceLocation",
            },
            "state": {"sequence": "G", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
    ),
    (
        "NM_181798.1:c.1007G>T",
        {
            "location": {
                "end": 1263,
                "sequenceReference": {
                    "refgetAccession": "SQ.KN07u-RFqd1dTyOWOG98HnOq87Nq-ZIg",
                    "type": "SequenceReference",
                },
                "start": 1262,
                "type": "SequenceLocation",
            },
            "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
    ),
    (
        "NC_000019.10:g.289464_289465insCACA",
        {
            "type": "Allele",
            "location": {
                "end": 289466,
                "sequenceReference": {
                    "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                    "type": "SequenceReference",
                },
                "start": 289464,
                "type": "SequenceLocation",
            },
            "state": {
                "length": 6,
                "repeatSubunitLength": 2,
                "sequence": "CACACA",
                "type": "ReferenceLengthExpression",
            },
        },
    ),
    (
        "NC_000019.10:g.289485_289500del",
        {
            "type": "Allele",
            "location": {
                "end": 289501,
                "sequenceReference": {
                    "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                    "type": "SequenceReference",
                },
                "start": 289480,
                "type": "SequenceLocation",
            },
            "state": {
                "length": 5,
                "repeatSubunitLength": 16,
                "sequence": "CGAGG",
                "type": "ReferenceLengthExpression",
            },
        },
    ),
)

hgvs_tests_to_hgvs_map = {
    "NC_000019.10:g.289464_289465insCACA": "NC_000019.10:g.289466_289467insCACA",
    "NC_000019.10:g.289485_289500del": "NC_000019.10:g.289486_289501del",
}


@pytest.mark.parametrize(("hgvsexpr", "expected"), hgvs_tests)
@pytest.mark.vcr
def test_hgvs(tlr, hgvsexpr, expected):
    # do_normalize defaults to true
    allele = tlr.translate_from(hgvsexpr, "hgvs")
    assert allele.model_dump(exclude_none=True) == expected

    to_hgvs = tlr.translate_to(allele, "hgvs")
    assert (hgvsexpr in to_hgvs) or (
        hgvs_tests_to_hgvs_map.get(hgvsexpr, hgvsexpr) in to_hgvs
    )


@pytest.mark.vcr
def test_rle_seq_limit(tlr):
    """Test that for an ReferenceLengthExpression over 50bp, the sequence is not
    returned when rle_seq_limit is set to 50bp, but is included in the state when
    rle_seq_limit is set to None.
    """
    # do_normalize defaults to true
    a_dict = {
        "location": {
            "end": 32331094,
            "sequenceReference": {
                "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                "type": "SequenceReference",
            },
            "start": 32331042,
            "type": "SequenceLocation",
        },
        "state": {
            "length": 104,
            "repeatSubunitLength": 52,
            "type": "ReferenceLengthExpression",
        },
        "type": "Allele",
    }
    input_hgvs_expr = "NC_000013.11:g.32331043_32331094dup"

    # use default rle_seq_limit
    allele_no_seq = tlr.translate_from(input_hgvs_expr, fmt="hgvs")
    assert allele_no_seq.model_dump(exclude_none=True) == a_dict

    with pytest.raises(
        AttributeError, match="'NoneType' object has no attribute 'root'"
    ):
        tlr.translate_to(allele_no_seq, "hgvs")

    # set rle_seq_limit to None
    allele_with_seq = tlr.translate_from(
        input_hgvs_expr, fmt="hgvs", rle_seq_limit=None
    )
    a_dict_with_seq = a_dict.copy()
    a_dict_with_seq["state"]["sequence"] = (
        "TTTAGTTGAACTACAGGTTTTTTTGTTGTTGTTGTTTTGATTTTTTTTTTTTTTTAGTTGAACTACAGGTTTTTTTGTTGTTGTTGTTTTGATTTTTTTTTTTT"
    )
    assert allele_with_seq.model_dump(exclude_none=True) == a_dict_with_seq

    output_hgvs_expr = tlr.translate_to(allele_with_seq, "hgvs")
    assert output_hgvs_expr == [input_hgvs_expr]


@pytest.mark.vcr
def test_to_hgvs_iri_ref_keyerror(tlr):
    # iriReference is passed
    iri_vo = models.Allele(
        **{  # noqa: PIE804
            "location": {
                "end": 1263,
                "start": 1262,
                "sequenceReference": "seqrefs.jsonc#/NM_181798.1",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        }
    )
    with pytest.raises(KeyError) as e:
        # even though the seqrefs.jsonc#/NM_181798.1 is a valid iri-reference for json schema, it is not yet handled here
        # we have to add functionality to address the handling of iri-references in the future
        tlr.translate_to(iri_vo, "hgvs")
    assert str(e.value) == "'ga4gh:seqrefs.jsonc#/NM_181798.1'"


@pytest.mark.vcr
def test_reference_allele_rle(tlr):
    """Test that reference alleles (REF==ALT) are normalized to ReferenceLengthExpression.

    Added to address https://github.com/ga4gh/vrs-python/issues/587
    """
    # Test with gnomad format
    gnomad_ref_allele = "1-100210778-AA-AA"
    allele = tlr._from_gnomad(gnomad_ref_allele)

    expected = {
        "type": "Allele",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
            },
            "start": 100210777,
            "end": 100210779,
        },
        "state": {
            "type": "ReferenceLengthExpression",
            "length": 2,
            "repeatSubunitLength": 2,
            "sequence": "AA",
        },
    }

    assert allele.model_dump(exclude_none=True) == expected

    # Test with SPDI format (REF==ALT)
    spdi_ref_allele = "NC_000001.11:100210777:AA:AA"
    allele_spdi = tlr._from_spdi(spdi_ref_allele)

    assert allele_spdi.model_dump(exclude_none=True) == expected

    # Test round-trip to SPDI
    to_spdi = tlr.translate_to(allele_spdi, "spdi", ref_seq_limit=None)
    assert len(to_spdi) == 1
    assert to_spdi[0] == spdi_ref_allele


# Microsatellite test cases for 21bp repeat unit
# https://github.com/ga4gh/vrs-python/discussions/592
# Tests deletion, insertion, and identity (no-change) variations
# Deletion/insertion: VOCA normalize to 930081-930152 with repeatSubunitLength=21
# Identity: Keep input coordinates 930089-930152 with repeatSubunitLength=63 (no VOCA normalization)
microsatellite_21bp_cases = [
    {
        "id": "deletion",
        "description": "Delete 1 copy from 3 copies (3->2 copies)",
        "hgvs": "NC_000001.11:g.930132_930152del",
        "spdi": "NC_000001.11:930081:GCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACC:GCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACC",
        "expected": {
            "type": "Allele",
            "location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                },
                "start": 930081,
                "end": 930152,
            },
            "state": {
                "type": "ReferenceLengthExpression",
                "length": 50,
                "sequence": "GCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACC",
                "repeatSubunitLength": 21,
            },
        },
    },
    {
        "id": "insertion",
        "description": "Insert 1 copy to 3 copies (3->4 copies)",
        "hgvs": "NC_000001.11:g.930152_930153insTTCCTCTCCTCCTGCCCCACC",
        "spdi": "NC_000001.11:930081:GCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACC:GCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACC",
        "expected": {
            "type": "Allele",
            "location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                },
                "start": 930081,
                "end": 930152,
            },
            "state": {
                "type": "ReferenceLengthExpression",
                "length": 92,
                "sequence": "GCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACC",
                "repeatSubunitLength": 21,
            },
        },
    },
    {
        "id": "identity",
        "description": "No change, 3 copies (same-as-ref does NOT do VOCA normalization)",
        "hgvs": "NC_000001.11:g.930090_930152=",
        "spdi": "NC_000001.11:930089:TTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACC:TTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACC",
        "expected": {
            "type": "Allele",
            "location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                },
                "start": 930089,
                "end": 930152,
            },
            "state": {
                "type": "ReferenceLengthExpression",
                "length": 63,
                "sequence": "TTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACCTTCCTCTCCTCCTGCCCCACC",
                "repeatSubunitLength": 63,
            },
        },
    },
]


@pytest.mark.parametrize("case", microsatellite_21bp_cases, ids=lambda c: c["id"])
@pytest.mark.vcr
def test_normalize_microsatellite_counts(tlr, case):
    """Test microsatellite deletion, insertion, and identity normalization behavior

    Tests three variations of a 21bp microsatellite:
    - Deletion and insertion: Apply VOCA normalization (roll left, find repeat unit)
    - Identity (same-as-ref): Do NOT apply VOCA normalization, use input coordinates

    https://github.com/ga4gh/vrs-python/discussions/592

    The microsatellite has a 21bp repeat unit when fully normalized.
    For deletion/insertion: normalize to 930081-930152 with repeatSubunitLength=21
    For identity: keep input coordinates (930089-930152) with repeatSubunitLength=63
    """
    # Test HGVS format
    allele_hgvs = tlr.translate_from(
        case["hgvs"], "hgvs", normalize=True, rle_seq_limit=100
    )
    assert allele_hgvs.model_dump(exclude_none=True) == case["expected"], (
        f"HGVS failed: {case['description']}"
    )

    # Test SPDI format
    allele_spdi = tlr.translate_from(
        case["spdi"], "spdi", normalize=True, rle_seq_limit=100
    )
    assert allele_spdi.model_dump(exclude_none=True) == case["expected"], (
        f"SPDI failed: {case['description']}"
    )


@pytest.mark.vcr
def test_translate_to_invalid_fmt(tlr):
    with pytest.raises(NotImplementedError, match="gnomad is not supported"):
        tlr.translate_to(models.Allele.model_validate(snv_output), fmt="gnomad")


# ---------------------------------------------------------------------------
# Regression tests for https://github.com/ga4gh/vrs-python/issues/364
# ("hgvs to vrs is returning valid results when hgvs has IncorrectReferenceAllele")
#
# These tests are hermetic: reference-sequence lookups are served from canned
# ground truth, so they run without seqrepo/UTA network access.


class _CannedDataProxy(_DataProxy):
    """Minimal data proxy backed by canned ground-truth reference sequences.

    Only the sequence/metadata lookups are faked; reference validation uses
    the real `_DataProxy.validate_ref_seq` logic.
    """

    # (accession, interbase start, interbase end) -> true reference sequence
    TRUTH: ClassVar[dict] = {
        # NM_006087.3 (TUBB4A): CDS is 373..1707, so c.900 == n.1272. The true
        # base there is G, as reported by the ClinGen Allele Registry for the
        # NM_006087.3:c.900C>A expression in issue #364 (IncorrectReferenceAllele:
        # "given=C, found=G"), independently confirmed against NCBI RefSeq.
        ("NM_006087.3", 1271, 1272): "G",
        # GRCh38 chr19:44908822, true ref C (matches the C>T test expression)
        ("NC_000019.10", 44908821, 44908822): "C",
    }

    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        return self.TRUTH[(identifier, start, end)]

    def get_metadata(self, _identifier: str) -> dict:
        return {"aliases": ["ga4gh:SQ." + "A" * 32], "length": 10**6}


def _c_to_n_nm006087(_self, sv):
    """Emulate the UTA c.->n. mapping for NM_006087.3.

    UTA is not reachable from every test environment, so apply the true
    mapping directly: the RefSeq CDS annotation for NM_006087.3 is 373..1707,
    hence c.900 maps to n.1272.
    """
    assert sv.ac == "NM_006087.3"
    offset = 372  # n. coordinate == c. coordinate + 372 on this transcript
    sv.posedit.pos.start.base += offset
    sv.posedit.pos.end.base += offset
    sv.type = "n"
    return sv


@pytest.fixture
def tlr_canned():
    """AlleleTranslator backed by canned reference data (no network/UTA)."""
    with (
        patch("hgvs.dataproviders.uta.connect", return_value=MagicMock()),
        patch.object(HgvsTools, "c_to_n", _c_to_n_nm006087),
    ):
        yield AlleleTranslator(data_proxy=_CannedDataProxy(), identify=False)


def test_from_hgvs_wrong_ref_allele_raises(tlr_canned):
    """A mismatched reference allele in the HGVS expression must raise, not
    silently produce a plausible-but-wrong VRS Allele (issue #364).
    """
    error_msg = (
        "Reference mismatch at NM_006087.3 position 1271-1272 "
        "(input gave 'C' but correct ref is 'G')"
    )

    with pytest.raises(DataProxyValidationError) as e:
        tlr_canned._from_hgvs("NM_006087.3:c.900C>A", do_normalize=False)
    assert str(e.value) == error_msg

    with pytest.raises(DataProxyValidationError) as e:
        tlr_canned.translate_from(
            "NM_006087.3:c.900C>A", fmt="hgvs", do_normalize=False
        )
    assert str(e.value) == error_msg


def test_from_hgvs_wrong_ref_allele_no_validation(tlr_canned):
    """require_validation=False keeps the legacy behavior: the allele is
    returned and the mismatch is only logged.
    """
    allele = tlr_canned._from_hgvs(
        "NM_006087.3:c.900C>A", do_normalize=False, require_validation=False
    )
    assert (allele.location.start, allele.location.end) == (1271, 1272)
    assert allele.state.sequence.root == "A"


def test_from_hgvs_correct_ref_allele_passes(tlr_canned):
    """HGVS expressions whose stated ref matches the reference sequence still
    translate cleanly (substitution and deletion-with-ref forms).
    """
    allele = tlr_canned.translate_from(
        "NC_000019.10:g.44908822C>T", fmt="hgvs", do_normalize=False
    )
    assert (allele.location.start, allele.location.end) == (44908821, 44908822)
    assert allele.state.sequence.root == "T"

    allele = tlr_canned.translate_from(
        "NC_000019.10:g.44908822delC", fmt="hgvs", do_normalize=False
    )
    assert (allele.location.start, allele.location.end) == (44908821, 44908822)
    assert allele.state.sequence.root == ""


def test_from_hgvs_wrong_ref_allele_del_raises(tlr_canned):
    """The deletion-with-ref form is validated too."""
    with pytest.raises(DataProxyValidationError) as e:
        tlr_canned.translate_from(
            "NC_000019.10:g.44908822delG", fmt="hgvs", do_normalize=False
        )
    assert "correct ref is 'C'" in str(e.value)


def test_from_hgvs_no_ref_allele_skips_validation(tlr_canned):
    """Edits that state no reference allele (ins/dup/bare del) skip validation
    entirely -- the canned proxy raises KeyError on any sequence lookup, so
    this fails if validation is attempted.
    """
    allele = tlr_canned.translate_from(
        "NC_000019.10:g.44908822_44908823insT", fmt="hgvs", do_normalize=False
    )
    assert (allele.location.start, allele.location.end) == (44908822, 44908822)
    assert allele.state.sequence.root == "T"
