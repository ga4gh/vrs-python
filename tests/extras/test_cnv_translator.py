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
    (
        "NC_000013.11:g.26440969_26443305del",
        models.CopyChange.EFO_0030069,
        {
            "copyChange": {"primaryCode": "EFO:0030069"},
            "digest": "bY_40M893gdWmNhU598V-T_dYtPJs-pp",
            "id": "ga4gh:CX.bY_40M893gdWmNhU598V-T_dYtPJs-pp",
            "location": {
                "digest": "4akcjXlbAu4xBKnxjOL_b4DM_20HOCA3",
                "end": 26443305,
                "id": "ga4gh:SL.4akcjXlbAu4xBKnxjOL_b4DM_20HOCA3",
                "sequenceReference": {
                    "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                    "type": "SequenceReference",
                },
                "start": 26440968,
                "type": "SequenceLocation",
            },
            "type": "CopyNumberChange",
        },
    ),
    (
        "NC_000013.11:g.32379315_32379819del",
        None,
        {
            "copyChange": {"primaryCode": "EFO:0030067"},
            "digest": "WKtHlbV6XCoKqvyeAJxdFj4ogw9ipDfQ",
            "id": "ga4gh:CX.WKtHlbV6XCoKqvyeAJxdFj4ogw9ipDfQ",
            "location": {
                "digest": "_TUGA9kX6JKdXzUklgN2zWkOvNu5pNmV",
                "end": 32379819,
                "id": "ga4gh:SL._TUGA9kX6JKdXzUklgN2zWkOvNu5pNmV",
                "sequenceReference": {
                    "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                    "type": "SequenceReference",
                },
                "start": 32379314,
                "type": "SequenceLocation",
            },
            "type": "CopyNumberChange",
        },
    ),
    (
        "NC_000013.11:g.32332787_32333388dup",
        models.CopyChange.EFO_0030071,
        {
            "copyChange": {"primaryCode": "EFO:0030071"},
            "digest": "EqI18-X9p8MUDr-Oz2J5GPQppEQKqWMU",
            "id": "ga4gh:CX.EqI18-X9p8MUDr-Oz2J5GPQppEQKqWMU",
            "location": {
                "digest": "UOA3zJOPfQxxRord_7pBkoMBpt46xcQq",
                "end": 32333388,
                "id": "ga4gh:SL.UOA3zJOPfQxxRord_7pBkoMBpt46xcQq",
                "sequenceReference": {
                    "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                    "type": "SequenceReference",
                },
                "start": 32332786,
                "type": "SequenceLocation",
            },
            "type": "CopyNumberChange",
        },
    ),
    (
        "NC_000013.11:g.32344743_32352093dup",
        None,
        {
            "copyChange": {"primaryCode": "EFO:0030070"},
            "digest": "eZVvwCSineeWTGKG3vbJvqkSUMW2JTCH",
            "id": "ga4gh:CX.eZVvwCSineeWTGKG3vbJvqkSUMW2JTCH",
            "location": {
                "digest": "17-a6p7m6QznwxZyVz2QdA9oPf5jTCyT",
                "end": 32352093,
                "id": "ga4gh:SL.17-a6p7m6QznwxZyVz2QdA9oPf5jTCyT",
                "sequenceReference": {
                    "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                    "type": "SequenceReference",
                },
                "start": 32344742,
                "type": "SequenceLocation",
            },
            "type": "CopyNumberChange",
        },
    ),
)


@pytest.mark.parametrize(("hgvsexpr" ,"copy_change", "expected"), from_hgvs_cx_tests)
@pytest.mark.vcr
def test_from_hgvs_cx(tlr, hgvsexpr, copy_change, expected):
    """Test that _from_hgvs works correctly for copy number change"""
    cx = tlr._from_hgvs(hgvsexpr, copy_change=copy_change)
    assert cx.model_dump(exclude_none=True) == expected


from_hgvs_cn_tests = (
    (
        "NC_000013.11:g.26440969_26443305del",
        1,
        {
            "copies": 1,
            "digest": "QnU97C-cRW431O9qWia9UCVRBDvGDH7I",
            "id": "ga4gh:CN.QnU97C-cRW431O9qWia9UCVRBDvGDH7I",
            "location": {
                "digest": "4akcjXlbAu4xBKnxjOL_b4DM_20HOCA3",
                "end": 26443305,
                "id": "ga4gh:SL.4akcjXlbAu4xBKnxjOL_b4DM_20HOCA3",
                "sequenceReference": {
                    "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                    "type": "SequenceReference",
                },
                "start": 26440968,
                "type": "SequenceLocation",
            },
            "type": "CopyNumberCount",
        },
    ),
    (
        "NC_000013.11:g.32332787_32333388dup",
        2,
        {
            "copies": 2,
            "digest": "a3aEKSzI46jK7_HRSpGNDap6Z1j_7kMM",
            "id": "ga4gh:CN.a3aEKSzI46jK7_HRSpGNDap6Z1j_7kMM",
            "location": {
                "digest": "UOA3zJOPfQxxRord_7pBkoMBpt46xcQq",
                "end": 32333388,
                "id": "ga4gh:SL.UOA3zJOPfQxxRord_7pBkoMBpt46xcQq",
                "sequenceReference": {
                    "refgetAccession": "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                    "type": "SequenceReference",
                },
                "start": 32332786,
                "type": "SequenceLocation",
            },
            "type": "CopyNumberCount",
        },
    ),
)


@pytest.mark.parametrize(("hgvsexpr", "copies" ,"expected"), from_hgvs_cn_tests)
@pytest.mark.vcr
def test_from_hgvs_cn(tlr, hgvsexpr, copies, expected):
    """Test that _from_hgvs works correctly for copy number count"""
    cn = tlr._from_hgvs(hgvsexpr, copies=copies)
    assert cn.model_dump(exclude_none=True) == expected
