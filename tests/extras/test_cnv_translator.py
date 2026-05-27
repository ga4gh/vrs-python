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
        models.CopyChange.COMPLETE_GENOMIC_LOSS,
        {
            "copyChange": "complete genomic loss",
            "digest": "eVlFPLW5eFChqk1GX5QBpHent0UchVhG",
            "id": "ga4gh:CX.eVlFPLW5eFChqk1GX5QBpHent0UchVhG",
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
            "copyChange": "loss",
            "digest": "j47pPHxYNH4pWr9tQQAE9uaQDgI8Lb2e",
            "id": "ga4gh:CX.j47pPHxYNH4pWr9tQQAE9uaQDgI8Lb2e",
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
        models.CopyChange.LOW_LEVEL_GAIN,
        {
            "copyChange": "low-level gain",
            "digest": "p0bDtEJrI7g1zIlQxlVprTwXyB9j2nGX",
            "id": "ga4gh:CX.p0bDtEJrI7g1zIlQxlVprTwXyB9j2nGX",
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
            "copyChange": "gain",
            "digest": "XkrTJHxkqmtHh3gDtkotuQoeJ0uwGSlL",
            "id": "ga4gh:CX.XkrTJHxkqmtHh3gDtkotuQoeJ0uwGSlL",
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


@pytest.mark.parametrize(("hgvsexpr", "copy_change", "expected"), from_hgvs_cx_tests)
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


@pytest.mark.parametrize(("hgvsexpr", "copies", "expected"), from_hgvs_cn_tests)
@pytest.mark.vcr
def test_from_hgvs_cn(tlr, hgvsexpr, copies, expected):
    """Test that _from_hgvs works correctly for copy number count"""
    cn = tlr._from_hgvs(hgvsexpr, copies=copies)
    assert cn.model_dump(exclude_none=True) == expected


@pytest.mark.vcr
def test_from_hgvs_cn_copies_zero(tlr):
    """Test that copies=0 produces CopyNumberCount, not CopyNumberChange.

    copies=0 is a valid input (homozygous deletion), but 0 is falsy in Python
    so it was previously treated as missing and fell through to CopyNumberChange.
    """
    cn = tlr._from_hgvs("NC_000013.11:g.26440969_26443305del", copies=0)
    assert cn.type == "CopyNumberCount"
    assert cn.copies == 0


@pytest.mark.vcr
def test_from_hgvs_coding_coordinates(tlr):
    """CDS-relative (``c.``) inputs are converted to transcript-relative
    coordinates via ``c_to_n`` before the VRS location is built; without the
    conversion the resulting ``SequenceLocation`` would point into the 5' UTR.
    """
    cx = tlr._from_hgvs("NM_001331029.1:c.100_200del")
    assert cx.model_dump(exclude_none=True) == {
        "id": "ga4gh:CX.lJiuo6QPsrSEKI7EeX0xpX5Iqx0R_Kal",
        "type": "CopyNumberChange",
        "digest": "lJiuo6QPsrSEKI7EeX0xpX5Iqx0R_Kal",
        "location": {
            "id": "ga4gh:SL.K0g9SmC4z2-ayJsTHFHuIAygGNv1UMgi",
            "type": "SequenceLocation",
            "digest": "K0g9SmC4z2-ayJsTHFHuIAygGNv1UMgi",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.MBIgVnoHFw34aFqNUVGM0zgjC3d-v8dK",
            },
            "start": 249,
            "end": 350,
        },
        "copyChange": "loss",
    }


# Uncertain-range HGVS expressions. See https://github.com/biocommons/hgvs/issues/225
# for the full list of ClinVar examples; this covers 4 balanced/one-sided/both-sided
# cases from that list. The remaining expressions (e.g. ClinVar 565301, 220591,
# 254064, 237630, 11692, 251062 with GRCh37 accessions) can be added here as more
# sequence data becomes available in the local seqrepo test fixture.
from_hgvs_uncertain_range_tests = (
    # Balanced uncertain range on both sides (ClinVar 11692)
    (
        "NC_000023.11:g.(133661675_133661730)_(133661850_133661926)del",
        {
            "copyChange": "loss",
            "digest": "5owdE9qzsgtoRCJjuPuyGOX1sf04yt_Q",
            "id": "ga4gh:CX.5owdE9qzsgtoRCJjuPuyGOX1sf04yt_Q",
            "location": {
                "digest": "PYFrFoP4BanWt5cE-VukSpVTtde4vCUM",
                "end": [133661850, 133661926],
                "id": "ga4gh:SL.PYFrFoP4BanWt5cE-VukSpVTtde4vCUM",
                "sequenceReference": {
                    "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                    "type": "SequenceReference",
                },
                "start": [133661674, 133661729],
                "type": "SequenceLocation",
            },
            "type": "CopyNumberChange",
        },
    ),
    # Left side unknown (ClinVar 425669)
    (
        "NC_000002.12:g.(?_202376935)_(202377551_202464808)del",
        {
            "copyChange": "loss",
            "digest": "BazsTr3EQi7sWmQcG0k9Pa2p9IMi-mDu",
            "id": "ga4gh:CX.BazsTr3EQi7sWmQcG0k9Pa2p9IMi-mDu",
            "location": {
                "digest": "fs4VrXLo-EoFBXN32RS7yMwHPaUeqpIz",
                "end": [202377551, 202464808],
                "id": "ga4gh:SL.fs4VrXLo-EoFBXN32RS7yMwHPaUeqpIz",
                "sequenceReference": {
                    "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                    "type": "SequenceReference",
                },
                "start": [None, 202376934],
                "type": "SequenceLocation",
            },
            "type": "CopyNumberChange",
        },
    ),
    # Right side unknown (ClinVar 425698)
    (
        "NC_000002.12:g.(202377551_202464808)_(202559947_?)del",
        {
            "copyChange": "loss",
            "digest": "X_DqftAyv19jqe5nFiUEdZamCwC7jwUl",
            "id": "ga4gh:CX.X_DqftAyv19jqe5nFiUEdZamCwC7jwUl",
            "location": {
                "digest": "mlA8Y2369Re_FK00IjFNwGRIFtRHAyLr",
                "end": [202559947, None],
                "id": "ga4gh:SL.mlA8Y2369Re_FK00IjFNwGRIFtRHAyLr",
                "sequenceReference": {
                    "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                    "type": "SequenceReference",
                },
                "start": [202377550, 202464807],
                "type": "SequenceLocation",
            },
            "type": "CopyNumberChange",
        },
    ),
    # Both sides unknown (ClinVar 220591)
    (
        "NC_000017.11:g.(?_58709859)_(58734342_?)del",
        {
            "copyChange": "loss",
            "digest": "KzAPXUkzFNk5KUOb8RyjdoPFZ8DCkfBt",
            "id": "ga4gh:CX.KzAPXUkzFNk5KUOb8RyjdoPFZ8DCkfBt",
            "location": {
                "digest": "WLMVDKSlxtyE9Py_t27VFaLrc8Uoqkme",
                "end": [58734342, None],
                "id": "ga4gh:SL.WLMVDKSlxtyE9Py_t27VFaLrc8Uoqkme",
                "sequenceReference": {
                    "refgetAccession": "SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
                    "type": "SequenceReference",
                },
                "start": [None, 58709858],
                "type": "SequenceLocation",
            },
            "type": "CopyNumberChange",
        },
    ),
)


@pytest.mark.parametrize(("hgvsexpr", "expected"), from_hgvs_uncertain_range_tests)
@pytest.mark.vcr
def test_from_hgvs_uncertain_range(tlr, hgvsexpr, expected):
    """Test uncertain-range HGVS expressions produce SequenceLocations with
    Range start/end (VRS 2.x).
    """
    cx = tlr._from_hgvs(hgvsexpr)
    assert cx.model_dump(exclude_none=True) == expected
