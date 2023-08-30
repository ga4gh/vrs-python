import pytest
from ga4gh.vrs import models, normalize

# >>> dp.get_sequence("refseq:NC_000019.10", 44908820, 44908830)
#  |820      |825      | 830
# ' G C G C C T G G C A '
#    |A| a1
#

allele_dict = {
    "location": {
        "end": 26090951,
        "start": 26090950,
        "sequence": "NC_000006.12",
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "C",
        "type": "LiteralSequenceExpression"
      },
    "type": "Allele"
}


allele_dict_sequence_reference = {
    "location": {
        "end": 26090951,
        "start": 26090950,
        "sequence": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV"
        },
        "type": "SequenceLocation"
    },
    "state": {
        "sequence": "C",
        "type": "LiteralSequenceExpression"
      },
    "type": "Allele"
}


allele_dict2 = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequence": "refseq:NC_000023.11",
        "start": [None, 155980375],
        "end": [155980377, None]
    },
    "state": {
        "sequence": "",
        "type": "LiteralSequenceExpression"
    }
}


allele_dict2_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequence": "refseq:NC_000023.11",
        "start": [None, 155980375],
        "end": [155980377, None]
    },
    "state": {
        "length": 0,
        "repeatSubunitLength": 2,
        "type": "ReferenceLengthExpression"
    }
}


allele_dict3 = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequence": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": [155980374, 155980375],
        "end": [155980377, 155980378]
    },
    "state": {
        "sequence": "",
        "type": "LiteralSequenceExpression"
    }
}


allele_dict4 = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequence": "refseq:NC_000023.11",
        "start": 155980373,
        "end": 155980375
    },
    "state": {
        "sequence": "GTGT",
        "type": "LiteralSequenceExpression"
    }
}

allele_dict4_normalized = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequence": "refseq:NC_000023.11",
        "start": 155980373,
        "end": 155980375
    },
    "state": {
        "length": 4,
        "repeatSubunitLength": 2,
        "sequence": "GTGT",
        "type": "ReferenceLengthExpression"
    }
}


@pytest.mark.vcr
def test_normalize_allele(rest_dataproxy):
    allele1 = models.Allele(**allele_dict)
    allele2 = normalize(allele1, rest_dataproxy)
    assert allele1 == allele2

    allele1_seq_ref = models.Allele(**allele_dict_sequence_reference)
    allele2_seq_ref = normalize(allele1_seq_ref, rest_dataproxy)
    assert allele1_seq_ref == allele2_seq_ref

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
