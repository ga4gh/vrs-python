import pytest


inputs = {
    "hgvs": "NC_000013.11:g.32936732G>C",
    "beacon": "13 : 32936732 G > C",
    "spdi": "NC_000013.11:32936731:1:C",
    "vcf": "13-32936732-G-C"
}

output = {
    'id': 'ga4gh:VA.n9ax-9x6gOC0OEt73VMYqCBfqfxG1XUH',
    'location': {
        'id': 'ga4gh:SL.v9K0mcjQVugxTDIcdi7GBJ_R6fZ1lsYq',
        'interval': {'end': 32936732, 'start': 32936731, 'type': 'SimpleInterval'},
        'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
        'type': 'SequenceLocation'},
    'state': {'sequence': 'C', 'type': 'SequenceState'},
    'type': 'Allele'
}



def test_from_beacon(tlr):
    assert tlr.from_beacon(inputs["beacon"]).as_dict() == output

def test_from_hgvs(tlr):
    assert tlr.from_hgvs(inputs["hgvs"]).as_dict() == output

def test_from_spdi(tlr):
    assert tlr.from_spdi(inputs["spdi"]).as_dict() == output

def test_from_vcf(tlr):
    assert tlr.from_vcf(inputs["vcf"]).as_dict() == output


    
hgvs_tests = (
    ("NC_000013.11:g.32936732C=", {'id': 'ga4gh:VA.n9ax-9x6gOC0OEt73VMYqCBfqfxG1XUH',
                                   'location': {'id': 'ga4gh:SL.v9K0mcjQVugxTDIcdi7GBJ_R6fZ1lsYq',
                                                'interval': {'end': 32936732,
                                                             'start': 32936731,
                                                             'type': 'SimpleInterval'},
                                                'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                                'type': 'SequenceLocation'},
                                   'state': {'sequence': 'C', 'type': 'SequenceState'},
                                   'type': 'Allele'}), 
    ("NC_000007.14:g.55181320A>T", {'id': 'ga4gh:VA.vU0meY5wGjpyRLCjSxCfb2Jlruyn2adL',
                                    'location': {'id': 'ga4gh:SL.5D9eG-ev4fA7mYIpOpDEe-4Am1lzPZlQ',
                                                 'interval': {'end': 55181320,
                                                              'start': 55181319,
                                                              'type': 'SimpleInterval'},
                                                 'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                                                 'type': 'SequenceLocation'},
                                    'state': {'sequence': 'T', 'type': 'SequenceState'},
                                    'type': 'Allele'}),
    ("NC_000007.14:g.55181220del", {'id': 'ga4gh:VA.csOXic4ezsVVEPJjM7jdcx4cCYuWNvFx',
                                    'location': {'id': 'ga4gh:SL.eDAO6enI-Mok9nCCJotVmsKzi0vwBF9t',
                                                 'interval': {'end': 55181220,
                                                              'start': 55181219,
                                                              'type': 'SimpleInterval'},
                                                 'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                                                 'type': 'SequenceLocation'},
                                    'state': {'sequence': '', 'type': 'SequenceState'},
                                    'type': 'Allele'}),
    ("NC_000007.14:g.55181230_55181231insGGCT", {'id': 'ga4gh:VA.mL71zVuJ7BKsB6U825nJuGv31S84puyd',
                                                 'location': {'id': 'ga4gh:SL.YRGVXC7g1ScsKl_z594KbS8FLflV3sLV',
                                                              'interval': {'end': 55181230,
                                                                           'start': 55181230,
                                                                           'type': 'SimpleInterval'},
                                                              'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                                                              'type': 'SequenceLocation'},
                                                 'state': {'sequence': 'GGCT', 'type': 'SequenceState'},
                                                 'type': 'Allele'}),
)
@pytest.mark.parametrize("hgvsexpr,expected", hgvs_tests)
def test_hgvs(tlr, hgvsexpr, expected):
    assert expected == tlr.from_hgvs(hgvsexpr).as_dict()
    
