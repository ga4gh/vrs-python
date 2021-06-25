from ga4gh.core import sha512t24u, ga4gh_digest, ga4gh_serialize, ga4gh_identify, is_pjs_instance
from ga4gh.vrs import models, vrs_deref, vrs_enref

allele_dict = {
    'location': {'interval': {
                     'end': 55181320,
                     'start': 55181319,
                     'type': 'SimpleInterval'},
                 'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                 'type': 'SequenceLocation'},
    'state': {'sequence': 'T', 'type': 'SequenceState'},
    'type': 'Allele'}

a = models.Allele(**allele_dict)


def test_vr(): 

    assert a.as_dict() == allele_dict

    assert is_pjs_instance(a.location)

    assert ga4gh_serialize(a.location.interval) == b'{"end":55181320,"start":55181319,"type":"SimpleInterval"}'

    assert ga4gh_serialize(a.location) == b'{"interval":{"end":55181320,"start":55181319,"type":"SimpleInterval"},"sequence_id":"F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceLocation"}'
    assert sha512t24u(ga4gh_serialize(a.location)) == '5D9eG-ev4fA7mYIpOpDEe-4Am1lzPZlQ'
    assert ga4gh_digest(a.location) == '5D9eG-ev4fA7mYIpOpDEe-4Am1lzPZlQ'
    assert ga4gh_identify(a.location) == 'ga4gh:VSL.5D9eG-ev4fA7mYIpOpDEe-4Am1lzPZlQ'

    assert ga4gh_serialize(a) == b'{"location":"5D9eG-ev4fA7mYIpOpDEe-4Am1lzPZlQ","state":{"sequence":"T","type":"SequenceState"},"type":"Allele"}'
    assert ga4gh_digest(a) == 'vU0meY5wGjpyRLCjSxCfb2Jlruyn2adL'
    assert ga4gh_identify(a) == 'ga4gh:VA.vU0meY5wGjpyRLCjSxCfb2Jlruyn2adL'

    assert a.as_dict() == {'location': {'interval': {'end': 55181320, 'start': 55181319, 'type': 'SimpleInterval'},
                                        'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                                        'type': 'SequenceLocation'},
                           'state': {'sequence': 'T', 'type': 'SequenceState'},
                           'type': 'Allele'}
    
    vros = {}
    a2 = vrs_enref(a, vros)
    assert ga4gh_identify(a) == ga4gh_identify(a2)
    assert a2.location == "ga4gh:VSL.5D9eG-ev4fA7mYIpOpDEe-4Am1lzPZlQ"
    assert a2.location in vros
    assert ga4gh_identify(a) in vros

    a3 = vrs_deref(a2, vros)
    assert a == a3
