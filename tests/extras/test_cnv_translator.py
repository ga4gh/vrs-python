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
     )
)


@pytest.mark.parametrize("hgvsexpr,copy_change,expected", from_hgvs_cx_tests)
@pytest.mark.vcr
def test_from_hgvs_cx(tlr, hgvsexpr ,copy_change, expected):
    """test that _from_hgvs works correctly for copy number change"""
    cx = tlr._from_hgvs(hgvsexpr, copy_change=copy_change)
    assert cx.model_dump(exclude_none=True) == expected


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
