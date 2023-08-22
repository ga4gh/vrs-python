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
        "sequence": "refseq:NC_000006.12",
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

seq_loc = {
        "type": "SequenceLocation",
        "sequence": "refseq:NC_000001.11",
        "start": [None, 244988599],
        "end": [244988601, None]
    }

allele_dict3 = {
    "type": "Allele",
    "location": seq_loc,
    "state": {
        "type": "RepeatedSequenceExpression",
        "seq_expr": {
            "location": seq_loc,
            "type": "DerivedSequenceExpression",
            "reverse_complement": False
        },
        "count": {
            "type": "Number",
            "value": 2
        }
    }
}


@pytest.mark.vcr
def test_normalize_allele(rest_dataproxy):
    allele1 = models.Allele(**allele_dict)
    allele2 = normalize(allele1, rest_dataproxy)
    assert allele1 == allele2

    allele1 = models.Allele(**allele_dict2)
    allele2 = normalize(allele1, rest_dataproxy)
    assert allele1 == allele2

    allele1 = models.Allele(**allele_dict3)
    allele2 = normalize(allele1, rest_dataproxy)
    assert allele1 == allele2
