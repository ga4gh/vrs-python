import pytest

inputs = {
    "hgvs": "NC_000013.11:g.32936732G>C",
    "beacon": "13 : 32936732 G > C",
    "spdi": "NC_000013.11:32936731:1:C",
    "gnomad": "13-32936732-G-C"
}

output = {
    'location': {
        'interval': {'end': 32936732, 'start': 32936731, 'type': 'SimpleInterval'},
        'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
        'type': 'SequenceLocation'},
    'state': {'sequence': 'C', 'type': 'SequenceState'},
    'type': 'Allele'
}



@pytest.mark.vcr
def test_from_beacon(tlr):
    assert tlr._from_beacon(inputs["beacon"]).as_dict() == output

@pytest.mark.vcr
def test_from_gnomad(tlr):
    assert tlr._from_gnomad(inputs["gnomad"]).as_dict() == output

@pytest.mark.vcr
def test_from_hgvs(tlr):
    assert tlr._from_hgvs(inputs["hgvs"]).as_dict() == output

@pytest.mark.vcr
def test_from_spdi(tlr):
    assert tlr._from_spdi(inputs["spdi"]).as_dict() == output



hgvs_tests = (
    ("NC_000013.11:g.32936732=",
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
@pytest.mark.vcr
def test_hgvs(tlr, hgvsexpr, expected):
    tlr.normalize = True
    allele = tlr.translate_from(hgvsexpr, "hgvs")
    assert expected == allele.as_dict()

    to_hgvs = tlr.translate_to(allele, "hgvs")
    assert 1 == len(to_hgvs)
    assert hgvsexpr == to_hgvs[0]


# TODO: Readd these tests
# @pytest.mark.vcr
# def test_errors(tlr):
#     with pytest.raises(ValueError):
#         tlr._from_beacon("bogus")
#
#     with pytest.raises(ValueError):
#         tlr._from_gnomad("NM_182763.2:c.688+403C>T")
#
#     with pytest.raises(ValueError):
#         tlr._from_hgvs("NM_182763.2:c.688+403C>T")
#
#     with pytest.raises(ValueError):
#         tlr._from_hgvs("NM_182763.2:c.688_690inv")
#
#     with pytest.raises(ValueError):
#         tlr._from_spdi("NM_182763.2:c.688+403C>T")

