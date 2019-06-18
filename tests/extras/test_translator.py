inputs = {
    "hgvs": "NC_000013.11:g.32936732G>C",
    "beacon": "13 : 32936732 G > C",
    "spdi": "NC_000013.11:32936731:1:C",
    "vcf": "13-32936732-G-C"
}

output = {
    'id': 'ga4gh:VA/cxH2DtdGp35-S0Eqj9Jf1-cH_yZiaf0U',
    'location': {
        'id': 'ga4gh:SL/HTlh4NTzBViUK0P5ZMinx-YJuNvlBmht',
        'interval': {'end': 32936732, 'start': 32936731, 'type': 'SimpleInterval'},
        'sequence_id': 'ga4gh:SQ/_0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
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
        
