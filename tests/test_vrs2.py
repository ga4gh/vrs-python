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
    assert a.model_dump(exclude_none=True) == allele_dict
    assert is_pydantic_instance(a)
    assert is_pydantic_instance(a.location)
    assert is_pydantic_instance(a.location.sequence)

    # Sequence Reference
    seqref = a.location.sequence
    seqref_serialized = ga4gh_serialize(seqref)
    assert seqref_serialized == b'{"refgetAccession":"SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceReference"}'
    assert sha512t24u(seqref_serialized) == 'OFEyBMeo55q3QRrxAY5FiDqnkdyf0GTV'
    assert ga4gh_digest(seqref) == 'OFEyBMeo55q3QRrxAY5FiDqnkdyf0GTV'
    assert ga4gh_identify(seqref) == 'ga4gh:SQR.OFEyBMeo55q3QRrxAY5FiDqnkdyf0GTV'

    # Location
    loc = a.location
    loc_serialized = ga4gh_serialize(loc)
    assert loc_serialized == b'{"end":55181320,"sequence":"OFEyBMeo55q3QRrxAY5FiDqnkdyf0GTV","start":55181319,"type":"SequenceLocation"}'
    assert sha512t24u(loc_serialized) == 'X0qrF7RfZxGIVIOddTYooZ_23D9mw6p6'
    assert ga4gh_digest(loc) == 'X0qrF7RfZxGIVIOddTYooZ_23D9mw6p6'
    assert ga4gh_identify(loc) == 'ga4gh:SL.X0qrF7RfZxGIVIOddTYooZ_23D9mw6p6'

    # Allele
    allele_serialized = ga4gh_serialize(a)
    assert allele_serialized == b'{"location":"X0qrF7RfZxGIVIOddTYooZ_23D9mw6p6","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}'
    assert sha512t24u(allele_serialized) == 'oFiLzDh37SoecjP7dceRaUfVlh32NnCg'
    assert ga4gh_digest(a) == 'oFiLzDh37SoecjP7dceRaUfVlh32NnCg'
    assert ga4gh_identify(a) == 'ga4gh:VA.oFiLzDh37SoecjP7dceRaUfVlh32NnCg'

    # Commenting out enref/deref tests.
    # We are deciding whether this will continue to be included in this library.
    # vros = {}
    # a2 = vrs_enref(a, vros)
    # assert ga4gh_identify(a) == ga4gh_identify(a2)
    # assert a2.location == "ga4gh:SL.Npx4j5beiNN9GSFTm8Ml6YxrNj_Ghkac"
    # assert a2.location in vros
    # assert ga4gh_identify(a) in vros

    # a3 = vrs_deref(a2, vros)
    # assert a == a3


def test_vr2_0():
    assert a.dict() == allele_dict


def test_iri():
    iri = models.IRI.model_construct("ga4gh:VA.asdf")
    assert is_curie_type(iri)
    assert iri.root == pydantic_copy(iri).root
    assert ga4gh_serialize(iri) == b'asdf'
