"""Out-of-bounds SequenceLocations are rejected on every translator input path

Sequence lengths used below:
    NM_000551.3     4560
    NP_001346993.1  193
    NC_000019.10    58617616
    NC_000007.14    159345973
    GRCh38:1        248956422
    NC_012920.1     16569
"""

import os
import re

import pytest

from ga4gh.vrs.dataproxy import DataProxyValidationError, SeqRepoRESTDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator, CnvTranslator

NC_000001_11 = "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
NC_000019_10 = "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
NM_000551_3 = "SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_"
NP_001346993_1 = "SQ.IPAWzkahAXVA3fBdoFluaU4NA3xTYUer"
NC_012920_1 = "SQ.k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct"


@pytest.fixture
def data_proxy() -> SeqRepoRESTDataProxy:
    # Metadata lookups are lru_cached per dataproxy instance, so a fresh instance per
    # test keeps each cassette self-contained regardless of test order
    return SeqRepoRESTDataProxy(
        base_url=os.environ.get("SEQREPO_REST_URL", "http://localhost:5000/seqrepo"),
        disable_healthcheck=True,
    )


@pytest.fixture
def allele_tlr(data_proxy: SeqRepoRESTDataProxy) -> AlleleTranslator:
    return AlleleTranslator(data_proxy=data_proxy)


@pytest.fixture
def cnv_tlr(data_proxy: SeqRepoRESTDataProxy) -> CnvTranslator:
    return CnvTranslator(data_proxy=data_proxy)


def _vrs_location(
    refget_accession: str,
    start: int | list[int | None],
    end: int | list[int | None],
) -> dict:
    return {
        "type": "SequenceLocation",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": refget_accession,
        },
        "start": start,
        "end": end,
    }


def _vrs_allele(
    refget_accession: str,
    start: int | list[int | None],
    end: int | list[int | None],
) -> dict:
    return {
        "type": "Allele",
        "location": _vrs_location(refget_accession, start, end),
        "state": {"type": "LiteralSequenceExpression", "sequence": "A"},
    }


def _bounds_msg(sequence_id: str, detail: str, seq_len: int) -> str:
    return (
        f"Location out of bounds on {sequence_id}: {detail} not within [0, {seq_len}]"
    )


OUT_OF_BOUNDS = [
    pytest.param(
        "allele_tlr",
        "hgvs",
        "NC_000019.10:g.58617617C>T",
        {},
        _bounds_msg(f"ga4gh:{NC_000019_10}", "end=58617617", 58617616),
        id="hgvs-g-past-end",
    ),
    pytest.param(
        "allele_tlr",
        "hgvs",
        "NM_000551.3:n.4561_4562insA",
        {},
        _bounds_msg(f"ga4gh:{NM_000551_3}", "start=4561, end=4561", 4560),
        id="hgvs-n-insertion-past-end",
    ),
    # ClinVar references the stop codon, which is not part of the protein sequence
    pytest.param(
        "allele_tlr",
        "hgvs",
        "NP_001346993.1:p.Ter194del",
        {},
        _bounds_msg(f"ga4gh:{NP_001346993_1}", "end=194", 193),
        id="hgvs-p-ter-at-length-plus-one",
    ),
    # Zero-width: an out-of-range fetch returns "" and would compare equal to the
    # empty reference, so only a coordinate check can catch this
    pytest.param(
        "allele_tlr",
        "spdi",
        "NM_000551.3:5000:0:AAA",
        {},
        _bounds_msg(f"ga4gh:{NM_000551_3}", "start=5000, end=5000", 4560),
        id="spdi-insertion-past-end",
    ),
    # Must report the bounds error, not "Reference mismatch ... correct ref is ''"
    pytest.param(
        "allele_tlr",
        "gnomad",
        "1-248956423-A-T",
        {},
        _bounds_msg("GRCh38:1", "end=248956423", 248956422),
        id="gnomad-past-end",
    ),
    pytest.param(
        "allele_tlr",
        "gnomad",
        "1-248956423-A-T",
        {"require_validation": False},
        _bounds_msg("GRCh38:1", "end=248956423", 248956422),
        id="gnomad-past-end-no-require-validation",
    ),
    pytest.param(
        "allele_tlr",
        "gnomad",
        "1-0-A-T",
        {},
        _bounds_msg("GRCh38:1", "start=-1", 248956422),
        id="gnomad-negative-start",
    ),
    pytest.param(
        "allele_tlr",
        "beacon",
        "1 : 248956423 A > T",
        {},
        _bounds_msg(f"ga4gh:{NC_000001_11}", "end=248956423", 248956422),
        id="beacon-past-end",
    ),
    pytest.param(
        "allele_tlr",
        "vrs",
        _vrs_allele(NM_000551_3, 99999999, 5),
        {},
        _bounds_msg(f"ga4gh:{NM_000551_3}", "start=99999999", 4560),
        id="vrs-allele-start-past-end-with-start-gt-end",
    ),
    pytest.param(
        "cnv_tlr",
        "hgvs",
        "NC_000007.14:g.159400000_159400100del",
        {},
        _bounds_msg("refseq:NC_000007.14", "start=159399999, end=159400100", 159345973),
        id="cnv-hgvs-copy-number-change-past-end",
    ),
]


IN_BOUNDS = [
    pytest.param(
        "allele_tlr",
        "hgvs",
        "NP_001346993.1:p.Leu193del",
        {"start": 192, "end": 193},
        id="hgvs-p-terminal-residue",
    ),
    pytest.param(
        "allele_tlr",
        "spdi",
        "NM_000551.3:4560:0:AAA",
        {"start": 4560, "end": 4560},
        id="spdi-insertion-at-end",
    ),
    pytest.param(
        "allele_tlr",
        "vrs",
        _vrs_allele(NC_012920_1, 16566, 5),
        {"start": 16566, "end": 5},
        id="vrs-allele-circular-start-gt-end",
    ),
]


@pytest.mark.parametrize(("tlr_fixture", "fmt", "var", "kwargs", "msg"), OUT_OF_BOUNDS)
@pytest.mark.vcr
def test_out_of_bounds(
    request: pytest.FixtureRequest,
    tlr_fixture: str,
    fmt: str,
    var: str | dict,
    kwargs: dict,
    msg: str,
) -> None:
    tlr = request.getfixturevalue(tlr_fixture)
    with pytest.raises(DataProxyValidationError, match=f"^{re.escape(msg)}$"):
        tlr.translate_from(var, fmt=fmt, **kwargs)


@pytest.mark.parametrize(("tlr_fixture", "fmt", "var", "expected_location"), IN_BOUNDS)
@pytest.mark.vcr
def test_in_bounds(
    request: pytest.FixtureRequest,
    tlr_fixture: str,
    fmt: str,
    var: str | dict,
    expected_location: dict,
) -> None:
    tlr = request.getfixturevalue(tlr_fixture)
    vo = tlr.translate_from(var, fmt=fmt)
    location = vo.location.model_dump()
    assert {k: location[k] for k in expected_location} == expected_location
