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


def test_normalize_other(dataproxy):
    i = models.SimpleInterval(start=5, end=6)
    assert i == normalize(i, dataproxy)
