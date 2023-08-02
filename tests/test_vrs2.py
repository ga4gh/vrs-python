from ga4gh.core import (
    sha512t24u,
    ga4gh_digest,
    ga4gh_serialize,
    ga4gh_identify,
    is_pydantic_instance,
    is_curie_type,
    pydantic_copy)
from ga4gh.vrs import models, vrs_deref, vrs_enref

allele_dict = {
    'location': {
        'end': 55181320,
        'start': 55181319,
        'sequence': {
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
    assert a.dict(exclude_none=True) == allele_dict
    assert is_pydantic_instance(a.location)

    assert ga4gh_serialize(
        a.location
    ) == b'{"end":55181320,"sequence":"OFEyBMeo55q3QRrxAY5FiDqnkdyf0GTV","start":55181319,"type":"SequenceLocation"}'
    # b'{"end":{"type":"Number","value":55181320},"sequence_id":"F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","start":{"type":"Number","value":55181319},"type":"SequenceLocation"}'
    assert sha512t24u(ga4gh_serialize(a.location)) == 'X0qrF7RfZxGIVIOddTYooZ_23D9mw6p6'
    assert ga4gh_digest(a.location) == 'X0qrF7RfZxGIVIOddTYooZ_23D9mw6p6'
    assert ga4gh_identify(a.location) == 'ga4gh:SL.X0qrF7RfZxGIVIOddTYooZ_23D9mw6p6'

    assert ga4gh_serialize(
        a
    ) == b'{"location":"Npx4j5beiNN9GSFTm8Ml6YxrNj_Ghkac","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}'
    assert ga4gh_digest(a) == 'RzhTjgnkCmLnaw3IBWnAubZs7eJHhho_'
    assert ga4gh_identify(a) == 'ga4gh:VA.RzhTjgnkCmLnaw3IBWnAubZs7eJHhho_'

    assert a.as_dict() == {
        'location': {
            'end': {'value': 55181320, 'type': 'Number'},
            'start': {'value': 55181319, 'type': 'Number'},
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
    assert a2.location == "ga4gh:SL.Npx4j5beiNN9GSFTm8Ml6YxrNj_Ghkac"
    assert a2.location in vros
    assert ga4gh_identify(a) in vros

    a3 = vrs_deref(a2, vros)
    assert a == a3


def test_vr2_0():
    assert a.dict() == allele_dict


def test_iri():
    iri = models.IRI.model_construct("ga4gh:VA.asdf")
    assert is_curie_type(iri)
    assert iri.root == pydantic_copy(iri).root
    assert ga4gh_serialize(iri) == b'asdf'
