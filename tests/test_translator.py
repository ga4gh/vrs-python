def test_from_beacon(tlr):
    assert tlr.from_beacon("13 : 32936732 G > C").as_dict() == {
        'location': {'interval': {'end': 32936732, 'start': 32936731,
        'type': 'SimpleInterval'}, 'sequence_id': 'GRCh38:13 ',
        'type': 'SequenceLocation'}, 'state': {'sequence': 'C',
        'type': 'SequenceState'}, 'type': 'Allele'}


def test_from_hgvs(tlr):
    assert tlr.from_hgvs("NM_012345.6:c.22A>T").as_dict() == {
        'location': { 'interval': {'end': 22, 'start': 21, 'type':
        'SimpleInterval'}, 'sequence_id': 'refseq:NM_012345.6',
        'type': 'SequenceLocation' }, 'state': {'sequence': 'T',
        'type': 'SequenceState'}, 'type': 'Allele' }


def test_from_spdi(tlr):
    assert tlr.from_spdi("NM_012345.6:21:1:T").as_dict() == {
          'location': { 'interval': {'end': 22, 'start': 21, 'type':
          'SimpleInterval'}, 'sequence_id': 'refseq:NM_012345.6',
          'type': 'SequenceLocation' }, 'state': {'sequence': 'T',
          'type': 'SequenceState'}, 'type': 'Allele' }


def test_from_vcf(tlr):
    assert tlr.from_vcf("1-55516888-G-GA").as_dict() == {'location':
        {'interval': {'end': 55516888, 'start': 55516887, 'type':
        'SimpleInterval'}, 'sequence_id': 'GRCh38:1', 'type':
        'SequenceLocation'}, 'state': {'sequence': 'GA', 'type':
        'SequenceState'}, 'type': 'Allele'}
        
