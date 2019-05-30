inputs = {
    "hgvs": "NC_000013.11:g.32936732G>C",
    "beacon": "13 : 32936732 G > C",
    "spdi": "NC_000013.11:32936731:1:C",
    "vcf": "13-32936732-G-C"
}

output = {
    'id': 'ga4gh:VAJBgSI1HBdpOYUNCWtRGwzhLtNrcdXAk8',
    'location': {
        'id': 'ga4gh:SL0FXQTd1CoM6ElQtD7qK1Ge6XGYhH6OZt',
        'interval': {'end': 32936732, 'start': 32936731, 'type': 'SimpleInterval'},
        'sequence_id': 'ga4gh:SQ_0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
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
        
