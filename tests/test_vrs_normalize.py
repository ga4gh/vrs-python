import pytest
from ga4gh.vrs import models, normalize

# >>> dp.get_sequence("refseq:NC_000019.10", 44908820, 44908830)
#  |820      |825      | 830
# ' G C G C C T G G C A '
#    |A| a1
#
allele_dict = {
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

# NC_000023.11:g.(?_155980375)_(155980377_?)del
uncertain_del_allele = {
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


# NC_000001.11:g.(?_244988599)_(244988601_?)dup
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

uncertain_dup_allele = {
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
def test_normalize_allele(dataproxy):
    allele1 = models.Allele(**allele_dict)
    allele2 = normalize(allele1, dataproxy)
    assert allele1 == allele2

    allele1.state.sequence = "C"
    allele3 = normalize(allele1, dataproxy)
    assert allele1 == allele3

    allele1.location.interval.end = 44908823
    allele1.state.sequence = ""
    allele4 = normalize(allele1, dataproxy)
    assert 44908820 == allele4.location.interval.start._value
    assert 44908824 == allele4.location.interval.end._value
    assert "GC" == allele4.state.sequence._value

    allele1 = models.Allele(**uncertain_del_allele)
    allele2 = normalize(allele1, dataproxy)
    assert allele1 == allele2

    allele1 = models.Allele(**uncertain_dup_allele)
    allele2 = normalize(allele1, dataproxy)
    assert allele1 == allele2


def test_normalize_other(dataproxy):
    i = models.SimpleInterval(start=5, end=6)
    assert i == normalize(i, dataproxy)
