from ga4gh.core import (
    sha512t24u,
    ga4gh_digest,
    ga4gh_serialize,
    ga4gh_identify,
    is_pydantic_instance,
    is_curie_type,
    pydantic_copy)
from ga4gh.vrs import models

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

allele_383650_dict = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequence": {
            "type": "SequenceReference",
            "refgetAccession": "ga4gh:SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI"
        },
        "start": 128325834,
        "end": 128325835
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "T"
    }
}
allele_417816_dict = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequence": {
            "type": "SequenceReference",
            "refgetAccession": "ga4gh:SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI"
        },
        "start": 128325809,
        "end": 128325810
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "T"
    }
}
allele_280320_dict = {
    "type": "Allele",
    "location": {
        "type": "SequenceLocation",
        "sequence": {
            "type": "SequenceReference",
            "refgetAccession": "ga4gh:SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI"
        },
        "start": 128322879,
        "end": 128322891
    },
    "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "G"
    }
}
allele_383650 = models.Allele(**allele_383650_dict)
allele_417816 = models.Allele(**allele_417816_dict)
allele_280320 = models.Allele(**allele_280320_dict)

haplotype_431012_dict = {
    "type": "Haplotype",
    "members": [allele_383650_dict, allele_417816_dict]
}
haplotype_431012 = models.Haplotype(**haplotype_431012_dict)

genotype_431013_dict = {
    "type": "Genotype",
    "count": 1,
    "members": [
        {
            "type": "GenotypeMember",
            "variation": haplotype_431012_dict,
            "count": 1
        },
        {
            "type": "GenotypeMember",
            "variation": allele_280320_dict,
            "count": 1
        }
    ]
}
genotype_431013 = models.Genotype(**genotype_431013_dict)


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


def test_haplotype():
    assert haplotype_431012.model_dump(exclude_none=True) == haplotype_431012_dict
    assert is_pydantic_instance(haplotype_431012)
    haplotype_serialized = ga4gh_serialize(haplotype_431012)
    assert haplotype_serialized == b'{"members":["WzQqNgeHveA5UCbiTWCivrVQEaSTKT5e","iHLnKCHfdlR6KB_OXLkzqCf86WGH91dq"],"type":"Haplotype"}'
    assert sha512t24u(haplotype_serialized) == 'i7MKrP2OU1nH0WQfwc4W_8DWVU0bac8X'
    assert ga4gh_digest(haplotype_431012) == 'i7MKrP2OU1nH0WQfwc4W_8DWVU0bac8X'
    assert ga4gh_identify(haplotype_431012) == 'ga4gh:HT.i7MKrP2OU1nH0WQfwc4W_8DWVU0bac8X'


def test_genotype():
    assert genotype_431013.model_dump(exclude_none=True) == genotype_431013_dict
    assert is_pydantic_instance(genotype_431013)
    genotype_serialized = ga4gh_serialize(genotype_431013)
    assert genotype_serialized == b'{"count":1,"members":[{"count":1,"type":"GenotypeMember","variation":"i7MKrP2OU1nH0WQfwc4W_8DWVU0bac8X"},{"count":1,"type":"GenotypeMember","variation":"DTN6tOnsmQ0aFodIDAzhkCzZUlfskiGb"}],"type":"Genotype"}'
    assert sha512t24u(genotype_serialized) == 'Wq0haS0V1I7fNlPHC5bpqU-tTr3xGREr'
    assert ga4gh_digest(genotype_431013) == 'Wq0haS0V1I7fNlPHC5bpqU-tTr3xGREr'
    assert ga4gh_identify(genotype_431013) == 'ga4gh:GT.Wq0haS0V1I7fNlPHC5bpqU-tTr3xGREr'


def test_iri():
    iri = models.IRI.model_construct("ga4gh:VA.asdf")
    assert is_curie_type(iri)
    assert iri.root == pydantic_copy(iri).root
    assert ga4gh_serialize(iri) == b'asdf'
