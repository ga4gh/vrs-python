import pytest

inputs = {
    "hgvs": "NC_000013.11:g.32936732G>C",
    "beacon": "13 : 32936732 G > C",
    "spdi": "NC_000013.11:32936731:1:C",
    "gnomad": "13-32936732-G-C"
}

output = {
    'location': {
        'interval': {
            'end': {'value': 32936732, 'type': 'Number'},
            'start': {'value': 32936731, 'type': 'Number'},
            'type': 'SequenceInterval'
        },
        'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
        'type': 'SequenceLocation'
    },
    'state': {
        'sequence': 'C',
        'type': 'LiteralSequenceExpression'
    },
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


@pytest.mark.vcr
def test_to_spdi(tlr):
    tlr.normalize = True
    spdiexpr = inputs["spdi"]
    allele = tlr.translate_from(spdiexpr, "spdi")
    to_spdi = tlr.translate_to(allele, "spdi")
    assert 1 == len(to_spdi)
    assert spdiexpr == to_spdi[0]

hgvs_tests = (
    ("NC_000013.11:g.32936732=", {
        'location': {
            'interval': {
                'end': {'value': 32936732, 'type': 'Number'},
                'start': {'value': 32936731, 'type': 'Number'},
                'type': 'SequenceInterval'
            },
            'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'C',
            'type': 'LiteralSequenceExpression'
        },
        'type': 'Allele'
    }),
    ("NC_000007.14:g.55181320A>T", {
        'location': {
            'interval': {
                'end': {'value': 55181320, 'type': 'Number'},
                'start': {'value': 55181319, 'type': 'Number'},
                'type': 'SequenceInterval'
            },
            'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'T',
            'type': 'LiteralSequenceExpression'
        },
        'type': 'Allele'
    }),
    ("NC_000007.14:g.55181220del", {
        'location': {
            'interval': {
                'end': {'value': 55181220, 'type': 'Number'},
                'start': {'value': 55181219, 'type': 'Number'},
                'type': 'SequenceInterval'
            },
            'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': '',
            'type': 'LiteralSequenceExpression'
        },
        'type': 'Allele'
    }),
    ("NC_000007.14:g.55181230_55181231insGGCT", {
        'location': {
            'interval': {
                'end': {'value': 55181230, 'type': 'Number'},
                'start': {'value': 55181230, 'type': 'Number'},
                'type': 'SequenceInterval'
            },
            'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'GGCT',
            'type': 'LiteralSequenceExpression'
        },
        'type': 'Allele'
    }),
    ("NC_000013.11:g.32331093_32331094dup", {
        'location': {
            'interval': {
                'end': {'value': 32331094, 'type': 'Number'},
                'start': {'value': 32331082, 'type': 'Number'},
                'type': 'SequenceInterval'
            },
            'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'TTTTTTTTTTTTTT',
            'type': 'LiteralSequenceExpression'
        },
        'type': 'Allele'
    }),
    ("NC_000013.11:g.32316467dup", {
        'location': {
            'interval': {
                'end': {'value': 32316467, 'type': 'Number'},
                'start': {'value': 32316466, 'type': 'Number'},
                'type': 'SequenceInterval'
            },
            'sequence_id': 'ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'AA',
            'type': 'LiteralSequenceExpression'
        },
        'type': 'Allele'
    }),
    ("NM_001331029.1:n.872A>G", {
        'location': {
            'interval': {
                'end': {'value': 872, 'type': 'Number'},
                'start': {'value': 871, 'type': 'Number'},
                'type': 'SequenceInterval'
            },
            'sequence_id': 'ga4gh:SQ.MBIgVnoHFw34aFqNUVGM0zgjC3d-v8dK',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'G',
            'type': 'LiteralSequenceExpression'
        },
        'type': 'Allele'
    }),
    ("NM_181798.1:n.1263G>T", {
        'location': {
            'interval': {
                'end': {'value': 1263, 'type': 'Number'},
                'start': {'value': 1262, 'type': 'Number'},
                'type': 'SequenceInterval'
            },
            'sequence_id': 'ga4gh:SQ.KN07u-RFqd1dTyOWOG98HnOq87Nq-ZIg',
            'type': 'SequenceLocation'
        },
        'state': {
            'sequence': 'T',
            'type': 'LiteralSequenceExpression'
        },
        'type': 'Allele'
    }),
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
