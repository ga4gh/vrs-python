import pytest

from vcr_support import vcr


inputs = {
    "hgvs": "NC_000013.11:g.32936732G>C",
    "beacon": "13 : 32936732 G > C",
    "spdi": "NC_000013.11:32936731:1:C",
    "vcf": "13-32936732-G-C"
}

output = {
    'location': {
        'interval': {'end': 32936732, 'start': 32936731, 'type': 'SimpleInterval'},
        'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
        'type': 'SequenceLocation'},
    'state': {'sequence': 'C', 'type': 'SequenceState'},
    'type': 'Allele'
}



@vcr.use_cassette()
def test_from_beacon(tlr):
    assert tlr.from_beacon(inputs["beacon"]).as_dict() == output

@vcr.use_cassette()
def test_from_hgvs(tlr):
    assert tlr.from_hgvs(inputs["hgvs"]).as_dict() == output

@vcr.use_cassette()
def test_from_spdi(tlr):
    assert tlr.from_spdi(inputs["spdi"]).as_dict() == output

@vcr.use_cassette()
def test_from_vcf(tlr):
    assert tlr.from_vcf(inputs["vcf"]).as_dict() == output


    
hgvs_tests = (
    ("NC_000013.11:g.32936732C=",
     {'location': {'interval': {'end': 32936732,
                                'start': 32936731,
                                'type': 'SimpleInterval'},
                   'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'C', 'type': 'SequenceState'},
      'type': 'Allele'}), 

    ("NC_000007.14:g.55181320A>T",
     {'location': {'interval': {'end': 55181320,
                                'start': 55181319,
                                'type': 'SimpleInterval'},
                   'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'T', 'type': 'SequenceState'},
      'type': 'Allele'}),

    ("NC_000007.14:g.55181220del",
     {'location': {'interval': {'end': 55181220,
                                'start': 55181219,
                                'type': 'SimpleInterval'},
                   'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                   'type': 'SequenceLocation'},
      'state': {'sequence': '', 'type': 'SequenceState'},
      'type': 'Allele'}),

    ("NC_000007.14:g.55181230_55181231insGGCT",
     {'location': {'interval': {'end': 55181230,
                                'start': 55181230,
                                'type': 'SimpleInterval'},
                   'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                   'type': 'SequenceLocation'},
      'state': {'sequence': 'GGCT', 'type': 'SequenceState'},
      'type': 'Allele'}),
)
@pytest.mark.parametrize("hgvsexpr,expected", hgvs_tests)
@vcr.use_cassette()
def test_hgvs(tlr, hgvsexpr, expected):
    assert expected == tlr.from_hgvs(hgvsexpr).as_dict()



def test_errors(tlr):
    with pytest.raises(ValueError):
        tlr.from_beacon("bogus")

    with pytest.raises(ValueError):
        tlr.from_hgvs("NM_182763.2:c.688+403C>T")
        
    with pytest.raises(ValueError):
        tlr.from_hgvs("NM_182763.2:c.688_690inv")
        
    with pytest.raises(ValueError):
        tlr.from_spdi("NM_182763.2:c.688+403C>T")
        
    with pytest.raises(ValueError):
        tlr.from_vcf("NM_182763.2:c.688+403C>T")
        
