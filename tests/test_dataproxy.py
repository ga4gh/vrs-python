import os
import re

import pytest

from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import (
    DataProxyValidationError,
    _DataProxy,
    create_dataproxy,
)


@pytest.mark.parametrize("dp", ["rest_dataproxy", "dataproxy"])
@pytest.mark.vcr
def test_data_proxies(dp, request):
    dataproxy = request.getfixturevalue(dp)
    r = dataproxy.get_metadata("NM_000551.3")
    assert r["length"] == 4560
    assert "ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_" in r["aliases"]

    r = dataproxy.get_metadata("NC_000013.11")
    assert r["length"] == 114364328
    assert "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT" in r["aliases"]

    seq = dataproxy.get_sequence("ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_")
    assert seq.startswith("CCTCGCCTCCGTTACAACGGCCTACGGTGCTGGAGGATCCTTCTGCGCAC")

    seq = dataproxy.get_sequence(
        "ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_", start=0, end=50
    )
    assert seq == "CCTCGCCTCCGTTACAACGGCCTACGGTGCTGGAGGATCCTTCTGCGCAC"


def test_data_proxy_configs():
    with pytest.raises(
        ValueError,
        match=re.escape(
            "create_dataproxy scheme must include provider (e.g., `seqrepo+http:...`)"
        ),
    ):
        create_dataproxy("file:///path/to/seqrepo/root")

    with pytest.raises(
        ValueError,
        match=re.escape("SeqRepo URI scheme seqrepo+fake-scheme not implemented"),
    ):
        create_dataproxy("seqrepo+fake-scheme://localhost:5000")

    with pytest.raises(
        ValueError, match="DataProxy provider fake-dataprovider not implemented"
    ):
        create_dataproxy("fake-dataprovider+http://localhost:5000")

    with pytest.raises(
        ValueError,
        match="No data proxy URI provided or found in GA4GH_VRS_DATAPROXY_URI",
    ):
        create_dataproxy(None)

    # check that fallback on env var works
    os.environ["GA4GH_VRS_DATAPROXY_URI"] = "seqrepo+:tests/data/seqrepo/latest"
    create_dataproxy(None)

    # check that method arg takes precedence by passing an invalid arg
    with pytest.raises(
        ValueError,
        match=re.escape(
            "create_dataproxy scheme must include provider (e.g., `seqrepo+http:...`)"
        ),
    ):
        create_dataproxy("file:///path/to/seqrepo/root")


class _StubDataProxy(_DataProxy):
    """Dataproxy serving only metadata, keyed by identifier"""

    def __init__(self, metadata: dict[str, dict]) -> None:
        self.metadata = metadata

    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        raise NotImplementedError

    def get_metadata(self, identifier: str) -> dict:
        return self.metadata[identifier]


BOUNDS_SEQ_ID = "refseq:NM_000551.3"
BOUNDS_REFGET_ID = "ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_"
BOUNDS_SEQ_LEN = 4560


def _bounds_dp(aliases: list[str]) -> _StubDataProxy:
    """Stub serving one sequence's metadata under each of its aliases"""
    md = {"length": BOUNDS_SEQ_LEN, "aliases": aliases}
    return _StubDataProxy(dict.fromkeys(aliases, md))


# Locations representable on a sequence of length BOUNDS_SEQ_LEN
location_bounds_valid_cases = [
    {"id": "interior", "start": 10, "end": 20},
    {"id": "zero-width-at-start", "start": 0, "end": 0},
    {"id": "terminal-residue", "start": 4559, "end": 4560},
    {"id": "insertion-at-end", "start": 4560, "end": 4560},
    {"id": "circular-start-gt-end", "start": 4000, "end": 5},
    {
        "id": "indefinite-end-open-upper",
        "start": 4400,
        "end": models.Range([4500, None]),
    },
    {
        "id": "indefinite-start-open-lower",
        "start": models.Range([None, 10]),
        "end": 20,
    },
    {
        "id": "definite-ranges",
        "start": models.Range([0, 10]),
        "end": models.Range([4500, 4560]),
    },
    {"id": "start-undefined", "start": None, "end": 10},
    {"id": "end-undefined", "start": 10, "end": None},
    {"id": "both-undefined", "start": None, "end": None},
]

# Locations not representable on a sequence of length BOUNDS_SEQ_LEN, with the
# offending coordinates as reported in the error message
location_bounds_invalid_cases = [
    {"id": "one-past-end", "start": 4559, "end": 4561, "detail": "end=4561"},
    {
        "id": "zero-width-past-end",
        "start": 5000,
        "end": 5000,
        "detail": "start=5000, end=5000",
    },
    {
        "id": "start-past-end-with-start-gt-end",
        "start": 99999999,
        "end": 5,
        "detail": "start=99999999",
    },
    {"id": "negative-start", "start": -1, "end": 1, "detail": "start=-1"},
    {"id": "negative-end", "start": 0, "end": -1, "detail": "end=-1"},
    {
        "id": "definite-end-past-end",
        "start": 4400,
        "end": models.Range([4500, 4600]),
        "detail": "end=[4500, 4600]",
    },
    {
        "id": "indefinite-end-lower-bound-past-end",
        "start": 4400,
        "end": models.Range([4561, None]),
        "detail": "end=[4561, None]",
    },
    {
        "id": "range-with-negative-member",
        "start": models.Range([-5, 10]),
        "end": 20,
        "detail": "start=[-5, 10]",
    },
]


@pytest.mark.parametrize("case", location_bounds_valid_cases, ids=lambda c: c["id"])
def test_validate_location_bounds_valid(case: dict) -> None:
    dp = _bounds_dp([BOUNDS_SEQ_ID, BOUNDS_REFGET_ID])
    dp.validate_location_bounds(BOUNDS_SEQ_ID, case["start"], case["end"])


@pytest.mark.parametrize("case", location_bounds_invalid_cases, ids=lambda c: c["id"])
def test_validate_location_bounds_invalid(case: dict) -> None:
    dp = _bounds_dp([BOUNDS_SEQ_ID, BOUNDS_REFGET_ID])
    expected_msg = (
        f"Location out of bounds on {BOUNDS_SEQ_ID} ({BOUNDS_REFGET_ID}): "
        f"{case['detail']} not within [0, {BOUNDS_SEQ_LEN}]"
    )
    with pytest.raises(DataProxyValidationError, match=f"^{re.escape(expected_msg)}$"):
        dp.validate_location_bounds(BOUNDS_SEQ_ID, case["start"], case["end"])


# Name used for the sequence in the error message, given the sequence_id passed in
# and the aliases of the sequence
location_bounds_sequence_name_cases = [
    {
        "id": "input-and-refget",
        "sequence_id": BOUNDS_SEQ_ID,
        "aliases": [BOUNDS_SEQ_ID, BOUNDS_REFGET_ID],
        "seq_name": f"{BOUNDS_SEQ_ID} ({BOUNDS_REFGET_ID})",
    },
    {
        "id": "bare-accession-coerced",
        "sequence_id": "NM_000551.3",
        "aliases": [BOUNDS_SEQ_ID, BOUNDS_REFGET_ID],
        "seq_name": f"{BOUNDS_SEQ_ID} ({BOUNDS_REFGET_ID})",
    },
    {
        "id": "refget-input-named-once",
        "sequence_id": BOUNDS_REFGET_ID,
        "aliases": [BOUNDS_SEQ_ID, BOUNDS_REFGET_ID],
        "seq_name": BOUNDS_REFGET_ID,
    },
    {
        "id": "no-refget-alias",
        "sequence_id": BOUNDS_SEQ_ID,
        "aliases": [BOUNDS_SEQ_ID],
        "seq_name": BOUNDS_SEQ_ID,
    },
]


@pytest.mark.parametrize(
    "case", location_bounds_sequence_name_cases, ids=lambda c: c["id"]
)
def test_validate_location_bounds_sequence_name(case: dict) -> None:
    dp = _bounds_dp(case["aliases"])
    expected_prefix = f"Location out of bounds on {case['seq_name']}: "
    with pytest.raises(
        DataProxyValidationError, match=f"^{re.escape(expected_prefix)}"
    ):
        dp.validate_location_bounds(case["sequence_id"], 0, BOUNDS_SEQ_LEN + 1)
