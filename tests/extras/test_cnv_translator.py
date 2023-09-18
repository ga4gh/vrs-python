import pytest

from ga4gh.vrs import models
from ga4gh.vrs.extras.translator import CnvTranslator


@pytest.fixture(scope="module")
def tlr(rest_dataproxy):
    return CnvTranslator(
        data_proxy=rest_dataproxy,
        default_assembly_name="GRCh38",
        identify=True,
        normalize=False,
    )



from_hgvs_cx_tests = (
    ("NC_000013.11:g.26440969_26443305del", models.CopyChange.EFO_0030069, {
        "id": "ga4gh:CX.UvIENZGNmCM5j38yHueYdd-BZLn-aLJM",
        "location": {
            "id": "ga4gh:SL.XBjSbazxe6h8lI14zfXIxqt6J718QEec",
            "start": 26440968,
            "end": 26443305,
            "sequenceReference": {
                "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "copyChange": "efo:0030069",
        "type": "CopyNumberChange"
    }),
    ("NC_000013.11:g.32379315_32379819del", None, {
        "id": "ga4gh:CX.uUOl9g6TKESU5pFKvHuXyqY8jUEhEVWj",
        "location": {
            "id": "ga4gh:SL.i_sproJkj7zUS3zrzzOtQH6QxjmWMZmz",
            "start": 32379314,
            "end": 32379819,
            "sequenceReference": {
                "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "copyChange": "efo:0030067",
        "type": "CopyNumberChange"
    }),
    ("NC_000013.11:g.32332787_32333388dup", models.CopyChange.EFO_0030071, {
        "id": "ga4gh:CX.lH7khAM3sm88wuXauCdIjGk5G-QXofwJ",
        "location": {
            "id": "ga4gh:SL.uXkfbPl7kyvoKNqN3Zy3LkeRP3LV0pEw",
            "start": 32332786,
            "end": 32333388,
            "sequenceReference": {
                "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "copyChange": "efo:0030071",
        "type": "CopyNumberChange"
    }),
    ("NC_000013.11:g.32344743_32352093dup", None, {
        "id": "ga4gh:CX.bDt5kf-fa40r7jycFslnA5wGYET_6s2J",
        "location": {
            "id": "ga4gh:SL.OsTsPEe27fx-GUtI1cFCi_DpAinjPApm",
            "start": 32344742,
            "end": 32352093,
            "sequenceReference": {
                "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "copyChange": "efo:0030070",
        "type": "CopyNumberChange"
    })
)


@pytest.mark.parametrize("hgvsexpr,copy_change,expected", from_hgvs_cx_tests)
@pytest.mark.vcr
def test_from_hgvs_cx(tlr, hgvsexpr ,copy_change, expected):
    """test that _from_hgvs works correctly for copy number change"""
    cx = tlr._from_hgvs(hgvsexpr, copy_change=copy_change)
    assert expected == cx.model_dump(exclude_none=True)


from_hgvs_cn_tests = (
    ("NC_000013.11:g.26440969_26443305del", 1, {
        "id": "ga4gh:CN.tPLSlUBW6j2dVHv0G81U_IFXt4Xawo7J",
        "location": {
            "id": "ga4gh:SL.XBjSbazxe6h8lI14zfXIxqt6J718QEec",
            "start": 26440968,
            "end": 26443305,
            "sequenceReference": {
                "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "copies": 1,
        "type": "CopyNumberCount"
    }),
    ("NC_000013.11:g.32332787_32333388dup", 2, {
        "id": "ga4gh:CN.OonsZ_O31GsaXwadjS759G9lbrzq_2SL",
        "location": {
            "id": "ga4gh:SL.uXkfbPl7kyvoKNqN3Zy3LkeRP3LV0pEw",
            "start": 32332786,
            "end": 32333388,
            "sequenceReference": {
                "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                "type": "SequenceReference"
            },
            "type": "SequenceLocation"
        },
        "copies": 2,
        "type": "CopyNumberCount"
    }),
)


@pytest.mark.parametrize("hgvsexpr,copies,expected", from_hgvs_cn_tests)
@pytest.mark.vcr
def test_from_hgvs_cn(tlr, hgvsexpr ,copies, expected):
    """test that _from_hgvs works correctly for copy number count"""
    cn = tlr._from_hgvs(hgvsexpr, copies=copies)
    assert expected == cn.model_dump(exclude_none=True)
