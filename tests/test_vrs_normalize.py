import pytest

from ga4gh.vrs import models, normalize

# >>> dp.get_sequence("refseq:NC_000019.10", 44908820, 44908830)
#  |820      |825      | 830
# ' G C G C C T G G C A '
#    |A| a1
#


allele_dict1 = {
    "location": {
        "end": 26090951,
        "start": 26090950,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV",
        },
        "type": "SequenceLocation",
    },
    "state": {"sequence": "C", "type": "LiteralSequenceExpression"},
    "type": "Allele",
}

allele_dict1_normalized = {
    "location": {
        "end": 26090951,
        "start": 26090950,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV",
        },
        "type": "SequenceLocation",
    },
    "state": {
        "type": "ReferenceLengthExpression",
        "length": 1,
        "sequence": "C",
        "repeatSubunitLength": 1,
    },
    "type": "Allele",
}


allele_dict2 = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [None, 155980375],
        "end": [155980377, None],
    },
    "state": {"sequence": "", "type": "LiteralSequenceExpression"},
}


allele_dict2_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [None, 155980375],
        "end": [155980377, None],
    },
    "state": {
        "length": 0,
        "repeatSubunitLength": 2,
        "type": "ReferenceLengthExpression",
    },
}


allele_dict3 = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [155980374, 155980375],
        "end": [155980377, 155980378],
    },
    "state": {"sequence": "", "type": "LiteralSequenceExpression"},
}


allele_dict4 = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": 155980373,
        "end": 155980375,
    },
    "state": {"sequence": "GTGT", "type": "LiteralSequenceExpression"},
}

allele_dict4_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": 155980373,
        "end": 155980375,
    },
    "state": {
        "length": 4,
        "repeatSubunitLength": 2,
        "sequence": "GTGT",
        "type": "ReferenceLengthExpression",
    },
}

allele_dict5 = {
    "location": {
        "end": 289464,
        "start": 289464,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
        },
        "type": "SequenceLocation",
    },
    "state": {"sequence": "CAGCAG", "type": "LiteralSequenceExpression"},
    "type": "Allele",
}

allele_dict5_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
        },
        "start": 289464,
        "end": 289469,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        "length": 11,
        "sequence": "CAGCAGCAGCA",
        "repeatSubunitLength": 3,
    },
}

# Another same-as-reference allele
allele_dict6 = {
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
        "type": "LiteralSequenceExpression",
        "sequence": "AA",
    },
}

allele_dict6_normalized = {
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


@pytest.mark.vcr
def test_normalize_allele(rest_dataproxy):
    # Test that a same-as-reference LSE is normalized to an RLE
    allele1 = models.Allele(**allele_dict1)
    allele2 = normalize(allele1, rest_dataproxy)
    assert allele2 == models.Allele(**allele_dict1_normalized)

    allele1 = models.Allele(**allele_dict2)
    allele2 = normalize(allele1, rest_dataproxy, rle_seq_limit=0)
    assert allele1 != allele2
    assert allele2 == models.Allele(**allele_dict2_normalized)

    # Definite ranges are not normalized
    allele3 = models.Allele(**allele_dict3)
    allele3_after_norm = normalize(allele3, rest_dataproxy)
    assert allele3_after_norm == allele3

    # Duplication
    allele4 = models.Allele(**allele_dict4)
    allele4_after_norm = normalize(allele4, rest_dataproxy)
    assert allele4_after_norm == models.Allele(**allele_dict4_normalized)

    # Duplication in non-integer-repeat ambiguous region
    allele5 = models.Allele(**allele_dict5)
    allele5_after_norm = normalize(allele5, rest_dataproxy)
    assert allele5_after_norm == models.Allele(**allele_dict5_normalized)

    # Same-as-reference allele (REF==ALT)
    # Added to address https://github.com/ga4gh/vrs-python/issues/587
    allele6 = models.Allele(**allele_dict6)
    allele6_after_norm = normalize(allele6, rest_dataproxy)
    assert allele6_after_norm == models.Allele(**allele_dict6_normalized)


# Simple deletion. ClinVar 3385321
# SPDI: NC_000001.11:66926:G:
clinvar_deletion = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 66926,
        "end": 66927,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "",
    },
}

clinvar_deletion_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 66926,
        "end": 66927,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        "length": 0,
        "repeatSubunitLength": 1,
    },
}

# Microsatellite deletion: ClinVar 4286633
# SPDI: NC_000001.11:766399:AATAAATA:AATA
clinvar_microsatellite = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 766400,
        "end": 766404,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "",
    },
}

clinvar_microsatellite_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 766399,
        "end": 766407,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        "length": 4,
        "repeatSubunitLength": 4,
    },
}

# Tandem repeat deletion. ClinVar 1658573
# SPDI: NC_000001.11:930136:CTCCTCCT:CTCCT
clinvar_tandem_repeat = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 930137,
        "end": 930140,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "",
    },
}

clinvar_tandem_repeat_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 930136,
        "end": 930144,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        "length": 5,
        "repeatSubunitLength": 3,
    },
}

# Repeat subunit insertion in microsatellite region. ClinVar 2672290
# SPDI: NC_000001.11:1752908:CCTCCTCCTCCTCCTCCTCCTCCTCCTC:CCTCCTCCTCCTCCTCCTCCTCCTCCTCCTC
# This tests the circular expansion logic for insertions in repeat regions
clinvar_microsatellite_insertion = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 1752908,
        "end": 1752908,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "CCT",
    },
}

clinvar_microsatellite_insertion_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 1752908,
        "end": 1752936,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        "length": 31,
        "repeatSubunitLength": 3,
    },
}


@pytest.mark.vcr
def test_normalize_clinvar_rle(rest_dataproxy):
    """Test normalization of ClinVar variants that should produce RLE.

    These test cases are pulled from ClinVar GRCh38 VCF.
    """
    # Simple deletion: AG>A (deletes G)
    deletion = models.Allele(**clinvar_deletion)
    deletion_norm = normalize(deletion, rest_dataproxy, rle_seq_limit=0)
    assert deletion_norm == models.Allele(**clinvar_deletion_normalized)

    # Microsatellite: GAATA>G (deletes AATA repeat)
    microsatellite = models.Allele(**clinvar_microsatellite)
    microsatellite_norm = normalize(microsatellite, rest_dataproxy, rle_seq_limit=0)
    assert microsatellite_norm == models.Allele(**clinvar_microsatellite_normalized)

    # Tandem repeat: TCTC>T (deletes CTC repeat)
    tandem_repeat = models.Allele(**clinvar_tandem_repeat)
    tandem_repeat_norm = normalize(tandem_repeat, rest_dataproxy, rle_seq_limit=0)
    assert tandem_repeat_norm == models.Allele(**clinvar_tandem_repeat_normalized)

    # Repeat subunit insertion in microsatellite region: C>CCCT (inserts CCT in CTC[10] region)
    microsatellite_insertion = models.Allele(**clinvar_microsatellite_insertion)
    microsatellite_insertion_norm = normalize(
        microsatellite_insertion, rest_dataproxy, rle_seq_limit=0
    )
    assert microsatellite_insertion_norm == models.Allele(
        **clinvar_microsatellite_insertion_normalized
    )
