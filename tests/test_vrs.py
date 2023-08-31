from ga4gh.core import sha512t24u, ga4gh_digest, ga4gh_serialize, ga4gh_identify, is_pydantic_instance
from ga4gh.vrs import models, vrs_deref, vrs_enref

allele_dict = {
    'location': {
        'end': 55181320,
        'start': 55181319,
        'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul'
        },
        # 'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
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

    # assert a.model_dump() == allele_dict # TODO with model_config['extra'] == allow this assertion will always fail

    assert is_pydantic_instance(a.location)

    assert ga4gh_serialize(
        a.location
    ) == b'{"end":55181320,"sequenceReference":"OFEyBMeo55q3QRrxAY5FiDqnkdyf0GTV","start":55181319,"type":"SequenceLocation"}'
    assert sha512t24u(ga4gh_serialize(a.location)) == 'TQ--CzByOOf6uc89QqNqEzrFLSCkNiv0'
    assert ga4gh_digest(a.location) == 'TQ--CzByOOf6uc89QqNqEzrFLSCkNiv0'
    assert ga4gh_identify(a.location) == 'ga4gh:SL.TQ--CzByOOf6uc89QqNqEzrFLSCkNiv0'

    assert ga4gh_serialize(
        a
    ) == b'{"location":"TQ--CzByOOf6uc89QqNqEzrFLSCkNiv0","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}'
    assert ga4gh_digest(a) == 'NqnE9Bl89TMfNYExrZbDINwbHlKP9_CT'
    assert ga4gh_identify(a) == 'ga4gh:VA.NqnE9Bl89TMfNYExrZbDINwbHlKP9_CT'

    # assert a.model_dump(exclude_none=True) == {
    #     'location': {
    #         'end': 55181320,
    #         'sequence': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
    #         'start': 55181319,
    #         'type': 'SequenceLocation'
    #     },
    #     'state': {
    #         'sequence': 'T',
    #         'type': 'LiteralSequenceExpression'
    #     },
    #     'type': 'Allele'
    # }

    # vros = {}
    # a2 = vrs_enref(a, vros)
    # assert ga4gh_identify(a) == ga4gh_identify(a2)
    # assert a2.location == "ga4gh:SL.Npx4j5beiNN9GSFTm8Ml6YxrNj_Ghkac"
    # assert a2.location in vros
    # assert ga4gh_identify(a) in vros
    #
    # a3 = vrs_deref(a2, vros)
    # assert a == a3
