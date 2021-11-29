from ga4gh.core import sha512t24u, ga4gh_digest, ga4gh_serialize, ga4gh_identify, is_pjs_instance
from ga4gh.vrs import models, vrs_deref, vrs_enref

allele_dict = {
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
}

a = models.Allele(**allele_dict)


def test_vr():

    assert a.as_dict() == allele_dict

    assert is_pjs_instance(a.location)

    assert ga4gh_serialize(a.location.interval) == b'{"end":{"type":"Number","value":55181320},"start":{"type":"Number","value":55181319},"type":"SequenceInterval"}'

    assert ga4gh_serialize(
        a.location
    ) == b'{"interval":{"end":{"type":"Number","value":55181320},"start":{"type":"Number","value":55181319},"type":"SequenceInterval"},"sequence_id":"F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceLocation"}'
    assert sha512t24u(ga4gh_serialize(a.location)) == '1tI-5In290Gtp1DOlMxO6z-K3u0OzlXl'
    assert ga4gh_digest(a.location) == '1tI-5In290Gtp1DOlMxO6z-K3u0OzlXl'
    assert ga4gh_identify(a.location) == 'ga4gh:VSL.1tI-5In290Gtp1DOlMxO6z-K3u0OzlXl'

    assert ga4gh_serialize(
        a
    ) == b'{"location":"1tI-5In290Gtp1DOlMxO6z-K3u0OzlXl","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}'
    assert ga4gh_digest(a) == '5Z7gWQGUuGAPe4Pw2_kJvnkhS2Q5jRhY'
    assert ga4gh_identify(a) == 'ga4gh:VA.5Z7gWQGUuGAPe4Pw2_kJvnkhS2Q5jRhY'

    assert a.as_dict() == {
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
    }

    vros = {}
    a2 = vrs_enref(a, vros)
    assert ga4gh_identify(a) == ga4gh_identify(a2)
    assert a2.location == "ga4gh:VSL.1tI-5In290Gtp1DOlMxO6z-K3u0OzlXl"
    assert a2.location in vros
    assert ga4gh_identify(a) in vros

    a3 = vrs_deref(a2, vros)
    assert a == a3
