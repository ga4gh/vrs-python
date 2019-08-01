from ga4gh.core import sha512t24u
from ga4gh.vr import models, ga4gh_digest, ga4gh_serialize, ga4gh_identify

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

    assert ga4gh_serialize(a.location.interval) == b'{"end":55181320,"start":55181319,"type":"SimpleInterval"}'

    assert ga4gh_serialize(a.location) == b'{"interval":{"end":55181320,"start":55181319,"type":"SimpleInterval"},"sequence_id":"F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceLocation"}'
    assert sha512t24u(ga4gh_serialize(a.location)) == '5D9eG-ev4fA7mYIpOpDEe-4Am1lzPZlQ'
    assert ga4gh_digest(a.location) == '5D9eG-ev4fA7mYIpOpDEe-4Am1lzPZlQ'
    assert ga4gh_identify(a.location) == 'ga4gh:VSL.5D9eG-ev4fA7mYIpOpDEe-4Am1lzPZlQ'

    assert ga4gh_serialize(a) == b'{"location":"5D9eG-ev4fA7mYIpOpDEe-4Am1lzPZlQ","state":{"sequence":"T","type":"SequenceState"},"type":"Allele"}'
    assert ga4gh_digest(a) == 'vU0meY5wGjpyRLCjSxCfb2Jlruyn2adL'
    assert ga4gh_identify(a) == 'ga4gh:VA.vU0meY5wGjpyRLCjSxCfb2Jlruyn2adL'

    assert a.as_dict() == {'_digest': 'vU0meY5wGjpyRLCjSxCfb2Jlruyn2adL',
                           'location': {'_digest': '5D9eG-ev4fA7mYIpOpDEe-4Am1lzPZlQ',
                                        'interval': {'end': 55181320, 'start': 55181319, 'type': 'SimpleInterval'},
                                        'sequence_id': 'ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul',
                                        'type': 'SequenceLocation'},
                           'state': {'sequence': 'T', 'type': 'SequenceState'},
                           'type': 'Allele'}
