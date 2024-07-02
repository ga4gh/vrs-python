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
      'digest': '5Zm6EIFzBHpoRpeMTKAiDlxEGmF_Cq3I',
      'id': 'ga4gh:CX.5Zm6EIFzBHpoRpeMTKAiDlxEGmF_Cq3I',
      'location': {'digest': 'CO7rlO4O2Dlsxw0KnNw4kmrJwwFXgOcj',
                   'end': 26443305,
                   'id': 'ga4gh:SL.CO7rlO4O2Dlsxw0KnNw4kmrJwwFXgOcj',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 26440968,
                   'type': 'SequenceLocation'},
      'type': 'CopyNumberChange'}),
    ("NC_000013.11:g.32379315_32379819del", None,
     {'copyChange': 'efo:0030067',
      'digest': '2LU_BbNdXeQf8sMAWaYOXX2_wTtywgr1',
      'id': 'ga4gh:CX.2LU_BbNdXeQf8sMAWaYOXX2_wTtywgr1',
      'location': {'digest': 'FX6u3sfgZSBh3H21T__fINBpZ3N6itPD',
                   'end': 32379819,
                   'id': 'ga4gh:SL.FX6u3sfgZSBh3H21T__fINBpZ3N6itPD',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 32379314,
                   'type': 'SequenceLocation'},
      'type': 'CopyNumberChange'}
     ),
    ("NC_000013.11:g.32332787_32333388dup", models.CopyChange.EFO_0030071,
     {'copyChange': 'efo:0030071',
      'digest': '0bOGlWb3yHF4Pl0JdrehXrwxhlIeGWyb',
      'id': 'ga4gh:CX.0bOGlWb3yHF4Pl0JdrehXrwxhlIeGWyb',
      'location': {'digest': 'Da2mvtG6UDi5DTIslWRZ0dzlxKqfRuQ7',
                   'end': 32333388,
                   'id': 'ga4gh:SL.Da2mvtG6UDi5DTIslWRZ0dzlxKqfRuQ7',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 32332786,
                   'type': 'SequenceLocation'},
      'type': 'CopyNumberChange'}
     ),
    ("NC_000013.11:g.32344743_32352093dup", None,
     {'copyChange': 'efo:0030070',
      'digest': '9ayJB8XFftenknAHA1yaPGB2v-J-bdNo',
      'id': 'ga4gh:CX.9ayJB8XFftenknAHA1yaPGB2v-J-bdNo',
      'location': {'digest': '8mY-AOOgME6kJxYgIKaMUPj6s1k55bMX',
                   'end': 32352093,
                   'id': 'ga4gh:SL.8mY-AOOgME6kJxYgIKaMUPj6s1k55bMX',
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
      'digest': '3wzMLF_0TdEadCyazYu39NtXu8UyGUsr',
      'id': 'ga4gh:CN.3wzMLF_0TdEadCyazYu39NtXu8UyGUsr',
      'location': {'digest': 'CO7rlO4O2Dlsxw0KnNw4kmrJwwFXgOcj',
                   'end': 26443305,
                   'id': 'ga4gh:SL.CO7rlO4O2Dlsxw0KnNw4kmrJwwFXgOcj',
                   'sequenceReference': {'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                                         'type': 'SequenceReference'},
                   'start': 26440968,
                   'type': 'SequenceLocation'},
      'type': 'CopyNumberCount'}
     ),
    ("NC_000013.11:g.32332787_32333388dup", 2,
     {'copies': 2,
      'digest': 'VK-0Kj1sZ0KNbkRVRW23N6FPjimIsUd9',
      'id': 'ga4gh:CN.VK-0Kj1sZ0KNbkRVRW23N6FPjimIsUd9',
      'location': {'digest': 'Da2mvtG6UDi5DTIslWRZ0dzlxKqfRuQ7',
                   'end': 32333388,
                   'id': 'ga4gh:SL.Da2mvtG6UDi5DTIslWRZ0dzlxKqfRuQ7',
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
