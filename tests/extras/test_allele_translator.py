import pytest
from pydantic import ValidationError

from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import DataProxyValidationError
from ga4gh.vrs.extras.translator import AlleleTranslator


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
        # two base expansion (9*A -> 11*A): subunit lengths 1 and 2 are both valid, and
        # VRS 2.1 selects the smallest (VRS 2.0 selected 2)
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
                "repeatSubunitLength": 1,
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
        allele_gnomad = tlr.to_vrs(gnomad, fmt="gnomad", rle_seq_limit=None)
        assert allele_gnomad.model_dump(exclude_none=True) == expected_dict
        allele_spdi = tlr.to_vrs(spdi, fmt="spdi", rle_seq_limit=None)
        assert allele_spdi.model_dump(exclude_none=True) == expected_dict

        allele_gnomad_to_spdi = tlr.from_vrs(
            allele_gnomad, fmt="spdi", ref_seq_limit=None
        )
        assert len(allele_gnomad_to_spdi) == 1
        assert allele_gnomad_to_spdi[0] == spdi


def test_from_invalid(tlr):
    with pytest.raises(
        ValueError, match="Unable to parse data as gnomad, hgvs, spdi, vrs"
    ):
        tlr.to_vrs("BRAF amplication")

    with pytest.raises(
        ValueError, match="Unable to parse data as gnomad, hgvs, spdi, vrs"
    ):
        tlr.to_vrs("BRAF amplication", assembly_name="GRCh37")


@pytest.mark.parametrize(
    ("expr", "fmt"),
    [
        # multi-allelic gnomAD-style input must not be truncated to its first ALT
        ("1-55516888-G-GA,GT", "gnomad"),
        ("19-44908822-C-T;foo", "gnomad"),
        ("NC_000019.10:44908821:1:T:extra", "spdi"),
    ],
)
def test_from_trailing_input_rejected(tlr, expr, fmt):
    """Expressions with trailing content must be rejected, not partially parsed"""
    with pytest.raises(ValueError, match=f"Unable to parse data as {fmt}"):
        tlr.to_vrs(expr, fmt=fmt)


@pytest.mark.vcr
def test_from_gnomad(tlr) -> None:
    """Public gnomAD translation preserves validation and normalization behavior."""
    do_normalize = False
    assert (
        tlr.from_gnomad(
            snv_inputs["gnomad"], do_normalize=do_normalize, unrelated_option=True
        ).model_dump(exclude_none=True)
        == snv_output
    )
    assert (
        tlr.from_gnomad(mito_inputs["gnomad"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == mito_output
    )
    assert (
        tlr.from_gnomad(
            deletion_inputs["gnomad"], do_normalize=do_normalize
        ).model_dump(exclude_none=True)
        == gnomad_deletion_output
    )
    assert (
        tlr.from_gnomad(
            insertion_inputs["gnomad"], do_normalize=do_normalize
        ).model_dump(exclude_none=True)
        == gnomad_insertion_output
    )
    assert (
        tlr.from_gnomad(
            duplication_inputs["gnomad"], do_normalize=do_normalize
        ).model_dump(exclude_none=True)
        == duplication_output
    )

    # do_normalize defaults to true
    assert (
        tlr.from_gnomad(snv_inputs["gnomad"]).model_dump(exclude_none=True)
        == snv_output
    )
    assert (
        tlr.from_gnomad(mito_inputs["gnomad"]).model_dump(exclude_none=True)
        == mito_output
    )
    assert (
        tlr.from_gnomad(deletion_inputs["gnomad"]).model_dump(exclude_none=True)
        == deletion_output_normalized
    )
    assert (
        tlr.from_gnomad(insertion_inputs["gnomad"]).model_dump(exclude_none=True)
        == insertion_output
    )
    assert (
        tlr.from_gnomad(duplication_inputs["gnomad"]).model_dump(exclude_none=True)
        == duplication_output_normalized
    )

    assert tlr.from_gnomad("17-83129587-GTTGWCACATGA-G")

    # Test valid characters
    assert tlr.from_gnomad(
        "7-2-ACGTURYKMSWBDHVN-ACGTURYKMSWBDHVN", require_validation=False
    )

    # Invalid input. Ref does not match regex
    assert not tlr.from_gnomad("13-32936732-helloworld-C")

    # Ref != Actual ref
    invalid_var = "13-32936732-G-C"
    error_msg = "Reference mismatch at GRCh38:13 position 32936731-32936732 (input gave 'G' but correct ref is 'C')"

    with pytest.raises(DataProxyValidationError) as e:
        tlr.from_gnomad(invalid_var)
    assert str(e.value) == error_msg

    with pytest.raises(DataProxyValidationError) as e:
        tlr.to_vrs(invalid_var, fmt="gnomad")
    assert str(e.value) == error_msg

    # require_validation set to False
    assert tlr.from_gnomad(invalid_var, require_validation=False)


@pytest.mark.vcr
def test_from_hgvs(tlr) -> None:
    """Translate HGVS expressions directly through the public method."""
    do_normalize = False
    assert (
        tlr.from_hgvs(
            snv_inputs["hgvs"], do_normalize=do_normalize, assembly_name="GRCh38"
        ).model_dump(exclude_none=True)
        == snv_output
    )
    assert (
        tlr.from_hgvs(mito_inputs["hgvs"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == mito_output
    )
    assert (
        tlr.from_hgvs(deletion_inputs["hgvs"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == deletion_output
    )
    assert (
        tlr.from_hgvs(insertion_inputs["hgvs"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == insertion_output
    )
    assert (
        tlr.from_hgvs(duplication_inputs["hgvs"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == duplication_output
    )


@pytest.mark.vcr
def test_from_spdi(tlr) -> None:
    """Translate SPDI expressions directly through the public method."""
    do_normalize = False
    assert (
        tlr.from_spdi(
            snv_inputs["spdi"], do_normalize=do_normalize, require_validation=False
        ).model_dump(exclude_none=True)
        == snv_output
    )
    assert (
        tlr.from_spdi(mito_inputs["spdi"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == mito_output
    )
    for spdi_del_expr in deletion_inputs["spdi"]:
        assert (
            tlr.from_spdi(spdi_del_expr, do_normalize=do_normalize).model_dump(
                exclude_none=True
            )
            == deletion_output
        ), spdi_del_expr
    for spdi_ins_expr in insertion_inputs["spdi"]:
        assert (
            tlr.from_spdi(spdi_ins_expr, do_normalize=do_normalize).model_dump(
                exclude_none=True
            )
            == insertion_output
        ), spdi_ins_expr
    assert (
        tlr.from_spdi(duplication_inputs["spdi"], do_normalize=do_normalize).model_dump(
            exclude_none=True
        )
        == duplication_output
    )


@pytest.mark.vcr
def test_to_spdi(tlr) -> None:
    """Direct public exports agree with format dispatch."""
    # do_normalize defaults to true
    spdiexpr = snv_inputs["spdi"]
    allele = tlr.to_vrs(spdiexpr, "spdi")
    to_spdi = tlr.from_vrs(allele, "spdi")
    assert tlr.to_spdi(allele, unrelated_option=True) == to_spdi
    assert len(to_spdi) == 1
    assert spdiexpr == to_spdi[0]


@pytest.mark.vcr
def test_to_spdi_with_ref(tlr) -> None:
    """Direct and dispatched SPDI exports honor the reference sequence limit."""
    spdi_expr_no_ref = "NC_000019.10:44908821:1:T"
    spdi_expr_with_ref = "NC_000019.10:44908821:C:T"

    allele_no_ref = tlr.to_vrs(spdi_expr_no_ref, "spdi")
    allele_with_ref = tlr.to_vrs(spdi_expr_with_ref, "spdi")
    assert allele_no_ref == allele_with_ref

    to_spdi_no_ref = tlr.from_vrs(allele_no_ref, "spdi", ref_seq_limit=0)
    assert len(to_spdi_no_ref) == 1
    assert spdi_expr_no_ref == to_spdi_no_ref[0]

    to_spdi_with_ref = tlr.from_vrs(allele_with_ref, "spdi", ref_seq_limit=None)
    assert tlr.to_spdi(allele_with_ref, ref_seq_limit=None) == to_spdi_with_ref
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
                "repeatSubunitLength": 1,
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
def test_hgvs(tlr, hgvsexpr, expected) -> None:
    """Direct public exports agree with format dispatch."""
    # do_normalize defaults to true
    allele = tlr.to_vrs(hgvsexpr, "hgvs")
    assert allele.model_dump(exclude_none=True) == expected

    to_hgvs = tlr.from_vrs(allele, "hgvs")
    assert tlr.to_hgvs(allele, ref_seq_limit=0) == to_hgvs
    assert (hgvsexpr in to_hgvs) or (
        hgvs_tests_to_hgvs_map.get(hgvsexpr, hgvsexpr) in to_hgvs
    )


@pytest.mark.vcr
def test_rle_seq_limit(tlr, monkeypatch) -> None:
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
    allele_no_seq = tlr.to_vrs(input_hgvs_expr, fmt="hgvs")
    assert allele_no_seq.model_dump(exclude_none=True) == a_dict

    with pytest.raises(
        AttributeError, match="'NoneType' object has no attribute 'root'"
    ):
        tlr.from_vrs(allele_no_seq, "hgvs")

    # set rle_seq_limit to None
    allele_with_seq = tlr.to_vrs(input_hgvs_expr, fmt="hgvs", rle_seq_limit=None)
    a_dict_with_seq = a_dict.copy()
    a_dict_with_seq["state"]["sequence"] = (
        "TTTAGTTGAACTACAGGTTTTTTTGTTGTTGTTGTTTTGATTTTTTTTTTTTTTTAGTTGAACTACAGGTTTTTTTGTTGTTGTTGTTTTGATTTTTTTTTTTT"
    )
    assert allele_with_seq.model_dump(exclude_none=True) == a_dict_with_seq

    # Omitted limits use the current instance setting; explicit None stays unlimited.
    monkeypatch.setattr(tlr, "rle_seq_limit", None)
    assert tlr.from_hgvs(input_hgvs_expr).state.sequence is not None
    assert tlr.from_hgvs(input_hgvs_expr, rle_seq_limit=0).state.sequence is None

    monkeypatch.setattr(tlr, "rle_seq_limit", 0)
    assert tlr.from_hgvs(input_hgvs_expr).state.sequence is None
    assert tlr.from_hgvs(input_hgvs_expr, rle_seq_limit=None).state.sequence is not None

    output_hgvs_expr = tlr.from_vrs(allele_with_seq, "hgvs")
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
        tlr.from_vrs(iri_vo, "hgvs")
    assert str(e.value) == "'ga4gh:seqrefs.jsonc#/NM_181798.1'"


@pytest.mark.vcr
def test_reference_allele_rle(tlr) -> None:
    """Test that reference alleles (REF==ALT) are normalized to ReferenceLengthExpression.

    Added to address https://github.com/ga4gh/vrs-python/issues/587
    """
    # Test with gnomad format
    gnomad_ref_allele = "1-100210778-AA-AA"
    allele = tlr.from_gnomad(gnomad_ref_allele)

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
    allele_spdi = tlr.from_spdi(spdi_ref_allele)

    assert allele_spdi.model_dump(exclude_none=True) == expected

    # Test round-trip to SPDI
    to_spdi = tlr.from_vrs(allele_spdi, "spdi", ref_seq_limit=None)
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
    allele_hgvs = tlr.to_vrs(case["hgvs"], "hgvs", normalize=True, rle_seq_limit=100)
    assert allele_hgvs.model_dump(exclude_none=True) == case["expected"], (
        f"HGVS failed: {case['description']}"
    )

    # Test SPDI format
    allele_spdi = tlr.to_vrs(case["spdi"], "spdi", normalize=True, rle_seq_limit=100)
    assert allele_spdi.model_dump(exclude_none=True) == case["expected"], (
        f"SPDI failed: {case['description']}"
    )


@pytest.mark.vcr
def test_translate_to_invalid_fmt(tlr):
    with pytest.raises(NotImplementedError, match="gnomad is not supported"):
        tlr.from_vrs(models.Allele.model_validate(snv_output), fmt="gnomad")


@pytest.mark.parametrize("direct", [False, True], ids=["dispatch", "direct"])
def test_from_vrs_dict(direct: bool) -> None:
    """Regression test for ga4gh/vrs-python#489.

    Translating a VRS dict must resolve the model class from its `type`
    instead of crashing with TypeError (the `models` module was subscripted).
    """
    tlr = AlleleTranslator(data_proxy=None, identify=False)

    # valid VRS dict translates to a VRS object (previously raised TypeError)
    allele = (
        tlr.from_vrs_dict(snv_output) if direct else tlr.to_vrs(snv_output, fmt="vrs")
    )
    assert isinstance(allele, models.Allele)
    assert allele.type == "Allele"
    assert allele.location.start == snv_output["location"]["start"]
    assert allele.location.end == snv_output["location"]["end"]

    # non-dict and missing-type inputs still return None
    assert tlr.from_vrs_dict("NC_000019.10:g.44908822C>T") is None
    assert tlr.from_vrs_dict({"location": {}}) is None


from_vrs_non_variation_type_cases = [
    {"id": "unknown-name", "type": "NotARealModel"},
    {"id": "non-variation-vrs-class", "type": "SequenceLocation"},
    {"id": "non-model-module-attribute", "type": "Field"},
    {"id": "variation-base-class", "type": "_VariationBase"},
    {"id": "non-string", "type": 5},
]


@pytest.mark.parametrize(
    "case", from_vrs_non_variation_type_cases, ids=lambda c: c["id"]
)
@pytest.mark.parametrize("direct", [False, True], ids=["dispatch", "direct"])
def test_from_vrs_non_variation_type(case: dict, direct: bool) -> None:
    """A VRS dict whose `type` is not a VRS variation class is rejected, rather than
    resolved to whatever `ga4gh.vrs.models` attribute has that name
    """
    tlr = AlleleTranslator(data_proxy=None, identify=False)
    parse = tlr.from_vrs_dict if direct else tlr.to_vrs
    options = {} if direct else {"fmt": "vrs"}

    with pytest.raises(ValidationError) as e:
        parse({"type": case["type"]}, **options)
    assert e.value.errors()[0]["type"] == "union_tag_invalid"


def test_from_hgvs_uncertain_range_raises(tlr) -> None:
    """Uncertain-range HGVS cannot be represented as an Allele with a literal
    sequence state; use CnvTranslator instead. See issue #609.

    No cassette needed: the rejection happens before any data-proxy call.
    """
    with pytest.raises(ValueError, match="Uncertain-range"):
        tlr.from_hgvs("NC_000019.9:g.(11211022_11213339)_(11217364_11218067)dup")


@pytest.mark.vcr
def test_to_hgvs_uncertain_range(tlr):
    """An Allele with a ``Range`` SequenceLocation.start/end round-trips to
    uncertain-range HGVS syntax. See issue #609.
    """
    allele = models.Allele(
        location=models.SequenceLocation(
            sequenceReference=models.SequenceReference(
                refgetAccession="SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g"
            ),
            start=models.Range(root=[None, 202376934]),
            end=models.Range(root=[202377551, 202464808]),
        ),
        state=models.LiteralSequenceExpression(sequence=models.sequenceString(root="")),
    )
    assert tlr.from_vrs(allele, "hgvs") == [
        "NC_000002.12:g.(?_202376935)_(202377551_202464808)del"
    ]
