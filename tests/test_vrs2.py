from pydantic import ValidationError
import pytest

from ga4gh.core import (
    sha512t24u,
    ga4gh_digest,
    ga4gh_serialize,
    ga4gh_identify,
    is_pydantic_instance,
    is_curie_type,
    pydantic_copy,
    use_ga4gh_compute_identifier_when,
    GA4GHComputeIdentifierWhen,
)
from ga4gh.vrs import models, vrs_enref, vrs_deref

allele_dict = {
    'location': {
        'end': 55181320,
        'start': 55181319,
        'sequenceReference': {
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
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI"
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
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI"
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
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI"
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
    assert is_pydantic_instance(a.location.sequenceReference)

    # Sequence Reference
    seqref = a.location.sequenceReference
    seqref_serialized = ga4gh_serialize(seqref)
    assert seqref_serialized is None
    assert ga4gh_digest(seqref) == 'F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul'
    assert ga4gh_identify(seqref) is None

    # Location
    loc = a.location
    loc_serialized = ga4gh_serialize(loc)
    assert loc_serialized == b'{"end":55181320,"sequenceReference":"F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","start":55181319,"type":"SequenceLocation"}'
    assert sha512t24u(loc_serialized) == 'wbBBFBTTyBcJPyjkK7z_dCcHFm5pE-2K'
    assert ga4gh_digest(loc) == 'wbBBFBTTyBcJPyjkK7z_dCcHFm5pE-2K'
    assert ga4gh_identify(loc) == 'ga4gh:SL.wbBBFBTTyBcJPyjkK7z_dCcHFm5pE-2K'

    # Allele
    allele_serialized = ga4gh_serialize(a)
    assert allele_serialized == b'{"location":"wbBBFBTTyBcJPyjkK7z_dCcHFm5pE-2K","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}'
    assert sha512t24u(allele_serialized) == 'hOZr7drvRxkUT_srSFVq1NCzvAJdKJlw'
    assert ga4gh_digest(a) == 'hOZr7drvRxkUT_srSFVq1NCzvAJdKJlw'
    assert ga4gh_identify(a) == 'ga4gh:VA.hOZr7drvRxkUT_srSFVq1NCzvAJdKJlw'

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

    with pytest.raises(ValidationError):
        models.Allele(**{
            "type": "Allele",
            "location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    # refgetAccession can't include a namespace prefix
                    "refgetAccession": "ga4gh:SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI"
                },
                "start": 128325834,
                "end": 128325835
            },
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": "T"
            }
        })
    with pytest.raises(ValidationError):
        models.Allele(**{
            "type": "Allele",
            "location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI"
                },
                "start": 128325834,
                "end": 128325835
            },
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": "T"
            },
            # digest can't include a namespace prefix
            "digest": "ga4gh:734G5mtNwe40do8F6GKuqQP4QxyjBqVp"
        })



def test_haplotype():
    assert haplotype_431012.model_dump(exclude_none=True) == haplotype_431012_dict
    assert is_pydantic_instance(haplotype_431012)
    haplotype_serialized = ga4gh_serialize(haplotype_431012)
    assert haplotype_serialized == b'{"members":["734G5mtNwe40do8F6GKuqQP4QxyjBqVp","bU3n0M2YQaV5C5ebODJYZ0GnbyCrOIHi"],"type":"Haplotype"}'
    assert sha512t24u(haplotype_serialized) == 'fFR5oRpeD8Cuq2hfs3bXd1rgJUQrQA26'
    assert ga4gh_digest(haplotype_431012) == 'fFR5oRpeD8Cuq2hfs3bXd1rgJUQrQA26'
    assert ga4gh_identify(haplotype_431012) == 'ga4gh:HT.fFR5oRpeD8Cuq2hfs3bXd1rgJUQrQA26'


def test_genotype():
    assert genotype_431013.model_dump(exclude_none=True) == genotype_431013_dict
    assert is_pydantic_instance(genotype_431013)
    genotype_serialized = ga4gh_serialize(genotype_431013)
    assert genotype_serialized == b'{"count":1,"members":[{"count":1,"type":"GenotypeMember","variation":"fFR5oRpeD8Cuq2hfs3bXd1rgJUQrQA26"},{"count":1,"type":"GenotypeMember","variation":"AUYSTKn2HElZ_Gg-Cv9Pm6Yx9Xpvx8Tm"}],"type":"Genotype"}'
    assert sha512t24u(genotype_serialized) == '51J0mMryCGjdce3qBpqNt4n_hXUQmw83'
    assert ga4gh_digest(genotype_431013) == '51J0mMryCGjdce3qBpqNt4n_hXUQmw83'
    assert ga4gh_identify(genotype_431013) == 'ga4gh:GT.51J0mMryCGjdce3qBpqNt4n_hXUQmw83'


def test_iri():
    iri = models.IRI.model_construct("ga4gh:VA.asdf")
    assert is_curie_type(iri)
    assert iri.root == pydantic_copy(iri).root
    assert ga4gh_serialize(iri) == b'asdf'


def test_enref():
    object_store = {}
    allele_383650_enreffed = vrs_enref(allele_383650, object_store=object_store)
    orig_no_loc = allele_383650.model_dump().copy()
    orig_no_loc.pop("location")
    actual_no_loc = allele_383650_enreffed.model_dump().copy()
    actual_no_loc.pop("location")
    assert orig_no_loc == actual_no_loc, "Original and enreffed match except for enreffed field"
    assert allele_383650_enreffed.location == 'ga4gh:SL.cWtFS2CsCI1E_ocNVu6PeFQaMtVxIE-L'
    assert allele_383650_enreffed.model_dump(exclude_none=True) == {
        'type': 'Allele',
        'location': 'ga4gh:SL.cWtFS2CsCI1E_ocNVu6PeFQaMtVxIE-L',
        'state': {
            'type': 'LiteralSequenceExpression',
            'sequence': 'T'}}


    dereffed = vrs_deref(allele_383650_enreffed, object_store=object_store)
    assert (dereffed.location.model_dump(exclude_none=True) == {
        'type': 'SequenceLocation',
        'sequenceReference': {
            'type': 'SequenceReference',
            'refgetAccession': 'SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI'
        },
        'start': 128325834,
        'end': 128325835})
    assert dereffed.location.model_dump(exclude_none=True) == allele_383650.location.model_dump(exclude_none=True)
    assert dereffed.model_dump() == allele_383650.model_dump()

def test_enref2():
    object_store = {}
    a = {
        "type": "Allele",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
            },
            "start": 44908821,
            "end": 44908822
        },
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "T"
        }
    }
    vo_a = models.Allele(**a)
    a_enreffed = vrs_enref(vo_a, object_store=object_store)
    orig_no_loc = vo_a.model_dump().copy()
    orig_no_loc.pop("location")
    actual_no_loc = a_enreffed.model_dump().copy()
    actual_no_loc.pop("location")
    assert orig_no_loc == actual_no_loc, "Original and enreffed match except for enreffed field"
    assert a_enreffed.location == 'ga4gh:SL.m4B7OEJ6J3q6gPakM8mSEqKZeQkj1KYC'
    assert a_enreffed.model_dump(exclude_none=True) == {
        'type': 'Allele',
        'location': 'ga4gh:SL.m4B7OEJ6J3q6gPakM8mSEqKZeQkj1KYC',
        'state': {
            'type': 'LiteralSequenceExpression',
            'sequence': 'T'
        }
    }

def test_class_refatt_map():
    class_refatt_map_expected = {
        'Allele': ['location'],
        'Haplotype': ['members'],
        '_CopyNumber': ['location'],
        'CopyNumberCount': ['location'],
        'CopyNumberChange': ['location'],
        'GenotypeMember': ['variation']
    }
    assert class_refatt_map_expected == models.class_refatt_map


def test_compute_identifiers_when():
    a = {
        "type": "Allele",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.jdEWLvLvT8827O59m1Agh5H3n6kTzBsJ",
            },
            "start": 44908821,
            "end": 44908822,
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "T"},
    }
    correct_id = "ga4gh:VA.PkeY9RbMt9CEFakQ0AgDdAQ7POUeoWR0"
    syntax_valid_id = "ga4gh:VA.39eae078d9bb30da2a5c5d1969cb1472"
    syntax_invalid_id = "ga4gh:12345"

    # when id property is missing
    vo_a = models.Allele(**a)
    assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.ALWAYS):
        assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.INVALID):
        assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.MISSING):
        assert ga4gh_identify(vo_a) == correct_id

    # when id property is none
    a["id"] = None
    vo_a = models.Allele(**a)
    assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.ALWAYS):
        assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.INVALID):
        assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.MISSING):
        assert ga4gh_identify(vo_a) == correct_id

    # when id property is blank
    a["id"] = ""
    vo_a = models.Allele(**a)
    assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.ALWAYS):
        assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.INVALID):
        assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.MISSING):
        assert ga4gh_identify(vo_a) == correct_id

    # when id property is syntactically invalid
    a["id"] = syntax_invalid_id
    vo_a = models.Allele(**a)
    assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.ALWAYS):
        assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.INVALID):
        assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.MISSING):
        assert ga4gh_identify(vo_a) == syntax_invalid_id

    # when id property is syntactically valid
    a["id"] = syntax_valid_id
    vo_a = models.Allele(**a)
    assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.ALWAYS):
        assert ga4gh_identify(vo_a) == correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.INVALID):
        assert ga4gh_identify(vo_a) == syntax_valid_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.MISSING):
        assert ga4gh_identify(vo_a) == syntax_valid_id

    # when id property is correct
    a["id"] = correct_id
    vo_a = models.Allele(**a)
    assert ga4gh_identify(vo_a) == correct_id
    assert ga4gh_identify(vo_a) is not correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.ALWAYS):
        assert ga4gh_identify(vo_a) == correct_id
        assert ga4gh_identify(vo_a) is not correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.INVALID):
        assert ga4gh_identify(vo_a) == correct_id
        assert ga4gh_identify(vo_a) is correct_id
    with use_ga4gh_compute_identifier_when(GA4GHComputeIdentifierWhen.MISSING):
        assert ga4gh_identify(vo_a) == correct_id
        assert ga4gh_identify(vo_a) is correct_id
