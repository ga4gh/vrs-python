import pytest

from ga4gh.vrs import models, normalize

# Single nucleotide same-as-reference allele.
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

# Ambiguous indefinite-outer 2 bp deletion. Should become RLE.
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

# Ambiguous definite 2-4bp deletion. Cannot be converted to RLE. (as opposed to allele_dict2 which can)
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

# Tandem duplication of GT, normalizes to RLE.
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

# Insertion of multiple repeat subunits ("CAG") into an existing repeating region.
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

# Another simple same-as-reference allele
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

# Multi-base substitution (step 2.b). ClinVar 1530016
# HGVS: NC_000001.11:g.939146_939147delinsTT
# VCF: 1-939146-GA-TT
# Substitutions remain as LiteralSequenceExpression (both ref and alt non-empty after trim)
clinvar_substitution_2bp = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 939145,
        "end": 939147,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "TT",
    },
}

clinvar_substitution_normalized_2bp = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 939145,
        "end": 939147,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "TT",
    },
}

# Unambiguous insertion in a repeat region (step 5.a).
# Insert "CGT" into a poly-A run at chr1:236900409-236900417.
# Terminal bases (C, T) don't match surrounding A's, so no rolling occurs.
# Should remain as LiteralSequenceExpression at the original position.
unambiguous_insertion_in_repeat = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 236900413,
        "end": 236900413,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "CGT",
    },
}

unambiguous_insertion_in_repeat_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 236900413,
        "end": 236900413,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "CGT",
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

    # Definite ambiguous ranges are not normalized
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

    # Multi-base substitution (step 2.b): both ref and alt non-empty after trim
    # Should remain as LiteralSequenceExpression, not converted to RLE
    substitution = models.Allele(**clinvar_substitution_2bp)
    substitution_norm = normalize(substitution, rest_dataproxy)
    assert substitution_norm == models.Allele(**clinvar_substitution_normalized_2bp)

    # Unambiguous insertion in a repeat region (step 5.a)
    # Terminal bases don't match context, so no rolling - stays at original position
    unambig_ins = models.Allele(**unambiguous_insertion_in_repeat)
    unambig_ins_norm = normalize(unambig_ins, rest_dataproxy)
    assert unambig_ins_norm == models.Allele(
        **unambiguous_insertion_in_repeat_normalized
    )


# Simple deletion from non-repeating region (no trim/rolling involved). ClinVar 3385321
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


###############################################################################
# Partial repeat insertion/deletion edge cases in CCT repeat region
# Region: chr1:1752908-1752936 = CCTCCTCCTCCTCCTCCTCCTCCTCCTC (9 full CCT units + trailing C)
# These tests verify behavior when insertions/deletions don't align with repeat boundaries
###############################################################################
#### MIDDLE INSERTIONS (around position 1752915, in unit 3) ####

# Insert "CT" (2 bases, < repeat unit) in middle of CCT region
# SPDI: NC_000001.11:1752915::CT
# Reference: CCT CCT C      CT CCT ...  (positions 1752908-...)
#                    ^-- insert "CT" at interbase position 1752915
# Variant:  CCT CCT C[CT]CT CCT ...
# The ref "CCT" (3bp) becomes "CCTCT" (5bp), which cycles with period 2 (CT CT C)
partial_repeat_insertion = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 1752915,
        "end": 1752915,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "CT",
    },
}

partial_repeat_insertion_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        # Expands to 3 ref bases (1752915-1752918: "CCT"), not the full CCT repeat region,
        # because the inserted "CT" creates a 2-base cycle pattern locally
        "start": 1752915,
        "end": 1752918,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        # ref "CCT" (3 bases) + ins "CT" (2 bases) = "CCTCT" (5 bases), cycles with period 2
        "length": 5,
        "repeatSubunitLength": 2,
    },
}


# Insert "CCTC" (4 bases, > repeat unit) in middle of CCT region
# SPDI: NC_000001.11:1752915::CCTC
# Reference: CCT CCT C      CT CCT CCT ...  (positions 1752908-...)
#                     ^-- insert "CCTC" at interbase position 1752915 (between C at char 1752914 and C at 1752915)
# Variant:   CCT CCT C[CCTC]CT CCT ...
# In the variant sequence, the subunit now identified is now CCTC (CCTC CCTC C)
# And the ref "CCTCC" (5bp) becomes "CCTCCCTCC" (9bp, period 4)
middle_ins_4bp = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 1752915,
        "end": 1752915,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "CCTC",
    },
}

middle_ins_4bp_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        # Rolls left to 1752911, expands to 5 ref bases (CCTCC)
        "start": 1752911,
        "end": 1752916,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        # ref "CCTCC" (5 bases) + ins "CCTC" (4 bases) = 9 bases, cycles with period 4
        "length": 9,
        "repeatSubunitLength": 4,
    },
}

#### TAIL INSERTIONS (at position 1752934, end of unit 9) ####

# Insert "CT" (2 bases, < repeat unit) at tail of CCT region
# SPDI: NC_000001.11:1752934::CT
# Reference: ...CCT CCT CCT C G A  (positions 1752926-1752938)
#                        ^-- insert "CT" at interbase position 1752934 (after T, before trailing C)
# Variant:   ...CCT CCT CCT[CT]C G A
# Result: Stays as LSE - at boundary where next base (C) doesn't continue the CT pattern
# No expansion occurs; output is the unchanged insertion
tail_ins_2bp = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 1752934,
        "end": 1752934,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "CT",
    },
}

tail_ins_2bp_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        # No expansion - stays at original position
        "start": 1752934,
        "end": 1752934,
    },
    "state": {
        # Stays as LSE because insertion at boundary doesn't form repeating pattern
        "type": "LiteralSequenceExpression",
        "sequence": "CT",
    },
}

# Insert "CCTC" (4 bases, > repeat unit) at tail of CCT region
# SPDI: NC_000001.11:1752934::CCTC
# Reference: ...CCT CCT CCT C G A  (positions 1752926-1752938)
#                        ^-- insert "CCTC" at interbase position 1752934 (after T, before trailing C)
# Variant:   ...CCT CCT CCT[CCTC]C G A → ...CCT CCT CCTCCCTC G A
# Left-aligns: ref "T" (1752933-1752934) + ins "CCTC" = "CCCTC" after normalization
# Result: Stays as LSE with sequence "CCCTC" - doesn't form repeating pattern at boundary
tail_ins_4bp = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 1752934,
        "end": 1752934,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "CCTC",
    },
}

tail_ins_4bp_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        # Rolls left by 1 to include the T at 1752933
        "start": 1752933,
        "end": 1752934,
    },
    "state": {
        # Stays as LSE - ref "T" + ins "CCTC" = "CCCTC" after left-alignment
        "type": "LiteralSequenceExpression",
        "sequence": "CCCTC",
    },
}

#### MIDDLE DELETIONS (around positions 1752912-1752916) ####

# Delete "CTCC" (4 bases, > repeat unit) in middle of CCT region
# SPDI: NC_000001.11:1752912:CTCC:
# Reference: CCT CCT CCT CCT CCT ...  (positions 1752908-...)
#               C[CTCC]T CCT ...      <-- delete positions 1752912-1752916
# Variant:   CCT C     T CCT ... → CCTCTCCT...
# Left-aligns to 1752911: ref span becomes "CCTCC" (5 bases)
# Ref "CCTCC" - del "CTCC" = 1 base remaining; CTCC has period 4, so repeatSubunitLength=4
deletion_spanning_boundary = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 1752912,
        "end": 1752916,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "",
    },
}

deletion_spanning_boundary_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        # Rolls left by 1 to 1752911; ref span becomes "CCTCC" (5 bases)
        # Doesn't expand to full CCT repeat region because deleted "CTCC" has period 4, not 3
        "start": 1752911,
        "end": 1752916,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        # ref "CCTCC" (5 bases) - del "CTCC" (4 bases) = 1 base remaining
        # repeatSubunitLength=4 reflects the deletion size, not the original CCT repeat
        "length": 1,
        "repeatSubunitLength": 4,
    },
}

# Delete "TC" (2 bases, < repeat unit) in middle of CCT region
# SPDI: NC_000001.11:1752913:TC:
# Reference: CCT CCT CCT CCT CCT ...  (positions 1752908-...)
#               T[TC]CT CCT ...       <-- delete positions 1752913-1752915
# Variant:   CCT T  CT CCT ... → CCTTCTCCT...
# Left-aligns to 1752912: ref span becomes "CTC" (3 bases)
# Ref "CTC" - del "TC" = 1 base remaining; TC has period 2, so repeatSubunitLength=2
middle_del_2bp = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 1752913,
        "end": 1752915,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "",
    },
}

middle_del_2bp_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        # Rolls left by 1 to 1752912; ref span becomes "CTC" (3 bases)
        "start": 1752912,
        "end": 1752915,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        # ref "CTC" (3 bases) - del "TC" (2 bases) = 1 base remaining
        "length": 1,
        "repeatSubunitLength": 2,
    },
}

#### TAIL DELETIONS (at positions 1752932-1752936) ####

# Delete "TC" (2 bases, < repeat unit) at tail of CCT region
# SPDI: NC_000001.11:1752934:TC:
# Reference: ...CCT CCT CCT C G A  (positions 1752926-1752938)
#                       T[TC]G A   <-- delete positions 1752934-1752936 (T from CCT + trailing C)
# Variant:   ...CCT CCT CC   G A → ...CCTCCTCCGA
# Left-aligns to 1752933: ref span becomes "TTC" (3 bases, positions 1752933-1752936)
# Ref "TTC" - del "TC" = 1 base remaining; TC has period 2, so repeatSubunitLength=2
tail_del_2bp = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 1752934,
        "end": 1752936,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "",
    },
}

tail_del_2bp_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        # Rolls left by 1 to 1752933; ref span becomes "TTC" (3 bases)
        "start": 1752933,
        "end": 1752936,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        # ref "TTC" (3 bases) - del "TC" (2 bases) = 1 base remaining
        "length": 1,
        "repeatSubunitLength": 2,
    },
}

# Delete "CCTC" (4 bases, > repeat unit) at tail of CCT region
# SPDI: NC_000001.11:1752932:CCTC:
# Reference: ...CCT CCT [CCTC] G A  (positions 1752926-1752938)
#                       ^-- delete positions 1752932-1752936 (CCTC = unit 9 CCT + trailing C)
# Variant:   ...CCT CCT       G A → ...CCTCCTGA
# Left-aligns: ref span stays at 1752932-1752936 = "CCTC" (4 bases)
# Ref "CCTC" - del "CCTC" = 0 bases remaining; CCTC has period 4, so repeatSubunitLength=4
tail_del_4bp = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        "start": 1752932,
        "end": 1752936,
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "",
    },
}

tail_del_4bp_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
        },
        # No rolling - stays at original position
        "start": 1752932,
        "end": 1752936,
    },
    "state": {
        "type": "ReferenceLengthExpression",
        # ref "CCTC" (4 bases) - del "CCTC" (4 bases) = 0 bases remaining
        "length": 0,
        "repeatSubunitLength": 4,
    },
}


@pytest.mark.vcr
def test_normalize_partial_rle_del_ins(rest_dataproxy):
    """Test normalization of partial repeat insertions/deletions in CCT repeat region.

    Tests insertions and deletions that don't align with the 3-base CCT repeat boundary.
    Region: chr1:1752908-1752936 = CCTCCTCCTCCTCCTCCTCCTCCTCCTC

    Cases tested:
    - Middle insertions: 2bp and 4bp insertions within the repeat region
    - Tail insertions: 2bp and 4bp insertions at the end of the repeat region
    - Middle deletions: 2bp and 4bp deletions within the repeat region
    - Tail deletions: 2bp and 4bp deletions at the end of the repeat region
    """
    # === MIDDLE INSERTIONS ===

    # Middle ins 2bp: Insert "CT" at 1752915 (already tested in partial_repeat_insertion)
    partial_ins = models.Allele(**partial_repeat_insertion)
    partial_ins_norm = normalize(partial_ins, rest_dataproxy, rle_seq_limit=0)
    assert partial_ins_norm == models.Allele(**partial_repeat_insertion_normalized)

    # Middle ins 4bp: Insert "CCTC" at 1752915
    # 4-base insertion into 3-base repeat creates period-4 RLE
    mid_ins_4 = models.Allele(**middle_ins_4bp)
    mid_ins_4_norm = normalize(mid_ins_4, rest_dataproxy, rle_seq_limit=0)
    assert mid_ins_4_norm == models.Allele(**middle_ins_4bp_normalized)

    # === TAIL INSERTIONS ===

    # Tail ins 2bp: Insert "CT" at 1752934
    # At boundary, stays as LSE (doesn't form repeating pattern)
    tail_ins_2 = models.Allele(**tail_ins_2bp)
    tail_ins_2_norm = normalize(tail_ins_2, rest_dataproxy, rle_seq_limit=0)
    assert tail_ins_2_norm == models.Allele(**tail_ins_2bp_normalized)

    # Tail ins 4bp: Insert "CCTC" at 1752934
    # At boundary, stays as LSE with left-aligned sequence
    tail_ins_4 = models.Allele(**tail_ins_4bp)
    tail_ins_4_norm = normalize(tail_ins_4, rest_dataproxy, rle_seq_limit=0)
    assert tail_ins_4_norm == models.Allele(**tail_ins_4bp_normalized)

    # === MIDDLE DELETIONS ===

    # Middle del 4bp: Delete "CTCC" at 1752912-1752916 (already tested in deletion_spanning_boundary)
    span_del = models.Allele(**deletion_spanning_boundary)
    span_del_norm = normalize(span_del, rest_dataproxy, rle_seq_limit=0)
    assert span_del_norm == models.Allele(**deletion_spanning_boundary_normalized)

    # Middle del 2bp: Delete "TC" at 1752913-1752915
    # 2-base deletion creates period-2 RLE
    mid_del_2 = models.Allele(**middle_del_2bp)
    mid_del_2_norm = normalize(mid_del_2, rest_dataproxy, rle_seq_limit=0)
    assert mid_del_2_norm == models.Allele(**middle_del_2bp_normalized)

    # === TAIL DELETIONS ===

    # Tail del 2bp: Delete "TC" at 1752934-1752936
    # At boundary, still becomes RLE with period 2
    tail_del_2 = models.Allele(**tail_del_2bp)
    tail_del_2_norm = normalize(tail_del_2, rest_dataproxy, rle_seq_limit=0)
    assert tail_del_2_norm == models.Allele(**tail_del_2bp_normalized)

    # Tail del 4bp: Delete "CCTC" at 1752932-1752936
    # Complete 4-base deletion, RLE with length=0
    tail_del_4 = models.Allele(**tail_del_4bp)
    tail_del_4_norm = normalize(tail_del_4, rest_dataproxy, rle_seq_limit=0)
    assert tail_del_4_norm == models.Allele(**tail_del_4bp_normalized)
