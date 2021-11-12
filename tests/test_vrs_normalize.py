import pytest
from ga4gh.vrs import models, normalize

# >>> dp.get_sequence("refseq:NC_000019.10", 44908820, 44908830)
#  |820      |825      | 830
# ' G C G C C T G G C A '
#    |A| a1
#
deprecated_allele_dict = {
    'location': {
        'interval': {
            'end': 44908822,
            'start': 44908821,
            'type': 'SimpleInterval'
        },
        'sequence_id': 'refseq:NC_000019.10',
        'type': 'SequenceLocation'
    },
    'state': {
        'sequence': 'A',
        'type': 'SequenceState'
    },
    'type': 'Allele'
}

allele_dict = {
    "location": {
        "interval": {
            "end": {
                "type": "Number",
                "value": 26090951
            },
            "start": {
                "type": "Number",
                "value": 26090950
            },
            "type": "SequenceInterval"
        },
        "sequence_id": "refseq:NC_000006.12",
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
        "sequence_id": "refseq:NC_000023.11",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "type": "IndefiniteRange",
                "comparator": "<=",
                "value": 155980375
            },
            "end": {
                "type": "IndefiniteRange",
                "comparator": ">=",
                "value": 155980377
            }
        }
    },
    "state": {
        "sequence": "",
        "type": "LiteralSequenceExpression"
    }
}

seq_loc = {
        "type": "SequenceLocation",
        "sequence_id": "refseq:NC_000001.11",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "type": "IndefiniteRange",
                "comparator": "<=",
                "value": 244988599
            },
            "end": {
                "type": "IndefiniteRange",
                "comparator": ">=",
                "value": 244988601
            }
        }
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
    allele1 = models.Allele(**deprecated_allele_dict)
    allele2 = normalize(allele1, rest_dataproxy)
    assert allele1 == allele2

    allele1.state.sequence = "C"
    allele3 = normalize(allele1, rest_dataproxy)
    assert allele1 == allele3

    allele1.location.interval.end = 44908823
    allele1.state.sequence = ""
    allele4 = normalize(allele1, rest_dataproxy)
    assert 44908820 == allele4.location.interval.start._value
    assert 44908824 == allele4.location.interval.end._value
    assert "GC" == allele4.state.sequence._value

    allele1 = models.Allele(**allele_dict)
    allele2 = normalize(allele1, rest_dataproxy)
    assert allele1 == allele2

    allele1 = models.Allele(**allele_dict2)
    allele2 = normalize(allele1, rest_dataproxy)
    assert allele1 == allele2

    allele1 = models.Allele(**allele_dict3)
    allele2 = normalize(allele1, rest_dataproxy)
    assert allele1 == allele2


def test_normalize_other(rest_dataproxy):
    i = models.SimpleInterval(start=5, end=6)
    assert i == normalize(i, rest_dataproxy)
