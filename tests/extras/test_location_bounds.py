"""Out-of-bounds SequenceLocations are rejected on every translator input path
except ``vrs``, whose input is trusted as-is

Sequence lengths used below:
    NM_000551.3     4560
    NP_001346993.1  193
    NC_000007.14    159345973
    GRCh38:1        248956422
"""

import os
import re

import pytest

from ga4gh.vrs.dataproxy import DataProxyValidationError, SeqRepoRESTDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator, CnvTranslator

# Refget accession each input sequence resolves to, named alongside it in errors
REFGET_ACCESSIONS = {
    "GRCh38:1": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
    "refseq:NC_000007.14": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
    "refseq:NM_000551.3": "ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_",
    "refseq:NP_001346993.1": "ga4gh:SQ.IPAWzkahAXVA3fBdoFluaU4NA3xTYUer",
}


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


def _bounds_msg(sequence_id: str, detail: str, seq_len: int) -> str:
    return (
        f"Location out of bounds on {sequence_id} ({REFGET_ACCESSIONS[sequence_id]}): "
        f"{detail} not within [0, {seq_len}]"
    )


out_of_bounds_cases = [
    {
        "id": "hgvs-n-insertion-past-end",
        "tlr_fixture": "allele_tlr",
        "fmt": "hgvs",
        "var": "NM_000551.3:n.4561_4562insA",
        "msg": _bounds_msg("refseq:NM_000551.3", "start=4561, end=4561", 4560),
    },
    # ClinVar references the stop codon, which is not part of the protein sequence
    {
        "id": "hgvs-p-ter-at-length-plus-one",
        "tlr_fixture": "allele_tlr",
        "fmt": "hgvs",
        "var": "NP_001346993.1:p.Ter194del",
        "msg": _bounds_msg("refseq:NP_001346993.1", "end=194", 193),
    },
    # Zero-width: an out-of-range fetch returns "" and would compare equal to the
    # empty reference, so only a coordinate check can catch this
    {
        "id": "spdi-insertion-past-end",
        "tlr_fixture": "allele_tlr",
        "fmt": "spdi",
        "var": "NM_000551.3:5000:0:AAA",
        "msg": _bounds_msg("refseq:NM_000551.3", "start=5000, end=5000", 4560),
    },
    # Must report the bounds error, not "Reference mismatch ... correct ref is ''"
    {
        "id": "gnomad-past-end",
        "tlr_fixture": "allele_tlr",
        "fmt": "gnomad",
        "var": "1-248956423-A-T",
        "msg": _bounds_msg("GRCh38:1", "end=248956423", 248956422),
    },
    {
        "id": "gnomad-past-end-no-require-validation",
        "tlr_fixture": "allele_tlr",
        "fmt": "gnomad",
        "var": "1-248956423-A-T",
        "kwargs": {"require_validation": False},
        "msg": _bounds_msg("GRCh38:1", "end=248956423", 248956422),
    },
    {
        "id": "gnomad-negative-start",
        "tlr_fixture": "allele_tlr",
        "fmt": "gnomad",
        "var": "1-0-A-T",
        "msg": _bounds_msg("GRCh38:1", "start=-1", 248956422),
    },
    {
        "id": "beacon-past-end",
        "tlr_fixture": "allele_tlr",
        "fmt": "beacon",
        "var": "1 : 248956423 A > T",
        "msg": _bounds_msg("GRCh38:1", "end=248956423", 248956422),
    },
    {
        "id": "cnv-hgvs-copy-number-change-past-end",
        "tlr_fixture": "cnv_tlr",
        "fmt": "hgvs",
        "var": "NC_000007.14:g.159400000_159400100del",
        "msg": _bounds_msg(
            "refseq:NC_000007.14", "start=159399999, end=159400100", 159345973
        ),
    },
]


in_bounds_cases = [
    {
        "id": "spdi-insertion-at-end",
        "tlr_fixture": "allele_tlr",
        "fmt": "spdi",
        "var": "NM_000551.3:4560:0:AAA",
        "expected_location": {"start": 4560, "end": 4560},
    },
]


@pytest.mark.parametrize("case", out_of_bounds_cases, ids=lambda c: c["id"])
@pytest.mark.vcr
def test_out_of_bounds(request: pytest.FixtureRequest, case: dict) -> None:
    tlr = request.getfixturevalue(case["tlr_fixture"])
    with pytest.raises(DataProxyValidationError, match=f"^{re.escape(case['msg'])}$"):
        tlr.translate_from(case["var"], fmt=case["fmt"], **case.get("kwargs", {}))


@pytest.mark.parametrize("case", in_bounds_cases, ids=lambda c: c["id"])
@pytest.mark.vcr
def test_in_bounds(request: pytest.FixtureRequest, case: dict) -> None:
    tlr = request.getfixturevalue(case["tlr_fixture"])
    location = tlr.translate_from(case["var"], fmt=case["fmt"]).location.model_dump()
    expected_location = case["expected_location"]
    assert {k: location[k] for k in expected_location} == expected_location
