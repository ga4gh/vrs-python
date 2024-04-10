import pytest

from ga4gh.vrs import models
from ga4gh.vrs.extras.translator import CnvTranslator


@pytest.fixture(scope="module")
def tlr(rest_dataproxy):
    return CnvTranslator(
        data_proxy=rest_dataproxy,
        default_assembly_name="GRCh38",
        identify=True,
    )



from_hgvs_cx_tests = (
    ("NC_000013.11:g.26440969_26443305del", models.CopyChange.EFO_0030069,
     {'copyChange': 'efo:0030069',
      'digest': 'TvAhuGK6HYf53mXoUnon60cZ7DC_UgM3',
      'id': 'ga4gh:CX.TvAhuGK6HYf53mXoUnon60cZ7DC_UgM3',
      'location': {'digest': '4akcjXlbAu4xBKnxjOL_b4DM_20HOCA3',
                   'end': 26443305,
                   'id': 'ga4gh:SL.4akcjXlbAu4xBKnxjOL_b4DM_20HOCA3',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 26440968,
                   'type': 'SequenceLocation'},
      'type': 'CopyNumberChange'}),
    ("NC_000013.11:g.32379315_32379819del", None,
     {'copyChange': 'efo:0030067',
      'digest': 'K1J2muiVvrDWqKnd5cMpFbTP0eJUfeuE',
      'id': 'ga4gh:CX.K1J2muiVvrDWqKnd5cMpFbTP0eJUfeuE',
      'location': {'digest': '_TUGA9kX6JKdXzUklgN2zWkOvNu5pNmV',
                   'end': 32379819,
                   'id': 'ga4gh:SL._TUGA9kX6JKdXzUklgN2zWkOvNu5pNmV',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 32379314,
                   'type': 'SequenceLocation'},
      'type': 'CopyNumberChange'}
     ),
    ("NC_000013.11:g.32332787_32333388dup", models.CopyChange.EFO_0030071,
     {'copyChange': 'efo:0030071',
      'digest': '5gODsVN83N1fe9Lc_Octy5rBlkYl8pGU',
      'id': 'ga4gh:CX.5gODsVN83N1fe9Lc_Octy5rBlkYl8pGU',
      'location': {'digest': 'UOA3zJOPfQxxRord_7pBkoMBpt46xcQq',
                   'end': 32333388,
                   'id': 'ga4gh:SL.UOA3zJOPfQxxRord_7pBkoMBpt46xcQq',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 32332786,
                   'type': 'SequenceLocation'},
      'type': 'CopyNumberChange'}
     ),
    ("NC_000013.11:g.32344743_32352093dup", None,
     {'copyChange': 'efo:0030070',
      'digest': '-sUe85R9UC_RxX7e_B1YsmsdyLsGvvmq',
      'id': 'ga4gh:CX.-sUe85R9UC_RxX7e_B1YsmsdyLsGvvmq',
      'location': {'digest': '17-a6p7m6QznwxZyVz2QdA9oPf5jTCyT',
                   'end': 32352093,
                   'id': 'ga4gh:SL.17-a6p7m6QznwxZyVz2QdA9oPf5jTCyT',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 32344742,
                   'type': 'SequenceLocation'},
      'type': 'CopyNumberChange'}
    ),
    (
        "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup",
        None,
        {
            "copyChange": "efo:0030070",
            "digest": "H0-_q06in6rsvLfq_5b-CSmP4ZQ6r7-Q",
            "id": "ga4gh:CX.H0-_q06in6rsvLfq_5b-CSmP4ZQ6r7-Q",
            "location": {
                "digest": "R3FeXqOiAu8Vms7QngINQwIxW904fdWY",
                "end": [33274278, 33417151],
                "id": "ga4gh:SL.R3FeXqOiAu8Vms7QngINQwIxW904fdWY",
                "sequenceReference": {
                    "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                    "type": "SequenceReference"
                },
                "start": [31060226, 31100350],
                "type": "SequenceLocation"
            },
            "type": "CopyNumberChange"
        }
    ),
    (
        "NC_000009.11:g.(?_108337304)_(108337428_?)del",
        None,
        {
            "copyChange": "efo:0030067",
            "digest": "ANthOqEGX8MIn0kXuyQcn9bouYfbFgjH",
            "id": "ga4gh:CX.ANthOqEGX8MIn0kXuyQcn9bouYfbFgjH",
            "location": {
                "digest": "lpDGeQvnz80iis8xSxoCX_Pulnu7wx4M",
                "end": [108337428, None],
                "id": "ga4gh:SL.lpDGeQvnz80iis8xSxoCX_Pulnu7wx4M",
                "sequenceReference": {
                    "refgetAccession": "SQ.HBckYGQ4wYG9APHLpjoQ9UUe9v7NxExt",
                    "type": "SequenceReference"
                },
                "start": [None, 108337303],
                "type": "SequenceLocation"
            },
            "type": "CopyNumberChange"
        }
    ),
    (
        "NC_000005.9:g.(90136803)_(90159675)dup",
        None,
        {
            "copyChange": "efo:0030070",
            "digest": "YcbXUe21Bt1wQDV7zGM0lacOupkxduFS",
            "id": "ga4gh:CX.YcbXUe21Bt1wQDV7zGM0lacOupkxduFS",
            "location": {
                "digest": "r82CARuf8IxOidMdvQCUcsXNp3XiHEVH",
                "end": 90159675,
                "id": "ga4gh:SL.r82CARuf8IxOidMdvQCUcsXNp3XiHEVH",
                "sequenceReference": {
                    "refgetAccession": "SQ.vbjOdMfHJvTjK_nqvFvpaSKhZillW0SX",
                    "type": "SequenceReference"
                },
                "start": 90136802,
                "type": "SequenceLocation"
            },
            "type": "CopyNumberChange"
        }
    ),
    (
        "NC_000009.11:g.108337304_(108337428_?)del",
        None,
        {
            "copyChange": "efo:0030067",
            "digest": "brfJaiKCnSw-mvc3K9sUIEAyCN620PuD",
            "id": "ga4gh:CX.brfJaiKCnSw-mvc3K9sUIEAyCN620PuD",
            "location": {
                "digest": "6myLdODZ8WgbEDXc3HLp88ZbG536NCM-",
                "end": [108337428, None],
                "id": "ga4gh:SL.6myLdODZ8WgbEDXc3HLp88ZbG536NCM-",
                "sequenceReference": {
                    "refgetAccession": "SQ.HBckYGQ4wYG9APHLpjoQ9UUe9v7NxExt",
                    "type": "SequenceReference"
                },
                "start": 108337303,
                "type": "SequenceLocation"
            },
            "type": "CopyNumberChange"
        }
    )
)


@pytest.mark.parametrize("hgvsexpr,copy_change,expected", from_hgvs_cx_tests)
@pytest.mark.vcr
def test_from_hgvs_cx(tlr, hgvsexpr ,copy_change, expected):
    """test that _from_hgvs works correctly for copy number change"""
    cx = tlr._from_hgvs(hgvsexpr, copy_change=copy_change)
    assert cx.model_dump(exclude_none=True) == expected

@pytest.mark.vcf
def test_from_hgvs_cx_invalid(tlr):
    """test that _from_hgvs works correctly for copy number change invalid input"""
    # Should fail since it's not g. or m.
    with pytest.raises(
        ValueError, match="Only 'g' and 'm' reference sequences are supported"
    ):
        tlr._from_hgvs("NM_001197320.1:c.281_283dup")


from_hgvs_cn_tests = (
    ("NC_000013.11:g.26440969_26443305del", 1,
     {'copies': 1,
      'digest': 'QnU97C-cRW431O9qWia9UCVRBDvGDH7I',
      'id': 'ga4gh:CN.QnU97C-cRW431O9qWia9UCVRBDvGDH7I',
      'location': {'digest': '4akcjXlbAu4xBKnxjOL_b4DM_20HOCA3',
                   'end': 26443305,
                   'id': 'ga4gh:SL.4akcjXlbAu4xBKnxjOL_b4DM_20HOCA3',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 26440968,
                   'type': 'SequenceLocation'},
      'type': 'CopyNumberCount'}
     ),
    ("NC_000013.11:g.32332787_32333388dup", 2,
     {'copies': 2,
      'digest': 'a3aEKSzI46jK7_HRSpGNDap6Z1j_7kMM',
      'id': 'ga4gh:CN.a3aEKSzI46jK7_HRSpGNDap6Z1j_7kMM',
      'location': {'digest': 'UOA3zJOPfQxxRord_7pBkoMBpt46xcQq',
                   'end': 32333388,
                   'id': 'ga4gh:SL.UOA3zJOPfQxxRord_7pBkoMBpt46xcQq',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 32332786,
                   'type': 'SequenceLocation'},
      'type': 'CopyNumberCount'}
     ),
)


@pytest.mark.parametrize("hgvsexpr,copies,expected", from_hgvs_cn_tests)
@pytest.mark.vcr
def test_from_hgvs_cn(tlr, hgvsexpr ,copies, expected):
    """test that _from_hgvs works correctly for copy number count"""
    cn = tlr._from_hgvs(hgvsexpr, copies=copies)
    assert cn.model_dump(exclude_none=True) == expected
