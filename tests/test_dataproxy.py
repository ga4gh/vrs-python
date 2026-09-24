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
    """Dataproxy serving only sequence lengths"""

    def __init__(self, lengths: dict[str, int]) -> None:
        self.lengths = lengths

    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        raise NotImplementedError

    def get_metadata(self, identifier: str) -> dict:
        return {"length": self.lengths[identifier], "aliases": []}


BOUNDS_SEQ_ID = "refseq:NM_000551.3"
BOUNDS_SEQ_LEN = 4560

# (start, end) that are representable on a sequence of length BOUNDS_SEQ_LEN
LOCATION_BOUNDS_VALID = [
    pytest.param(10, 20, id="interior"),
    pytest.param(0, 0, id="zero-width-at-start"),
    pytest.param(4559, 4560, id="terminal-residue"),
    pytest.param(4560, 4560, id="insertion-at-end"),
    pytest.param(4000, 5, id="circular-start-gt-end"),
    pytest.param(4400, models.Range([4500, None]), id="indefinite-end-open-upper"),
    pytest.param(models.Range([None, 10]), 20, id="indefinite-start-open-lower"),
    pytest.param(
        models.Range([0, 10]), models.Range([4500, 4560]), id="definite-ranges"
    ),
    pytest.param(None, 10, id="start-undefined"),
    pytest.param(10, None, id="end-undefined"),
    pytest.param(None, None, id="both-undefined"),
]

# (start, end, offending coordinates as reported in the error message)
LOCATION_BOUNDS_INVALID = [
    pytest.param(4559, 4561, "end=4561", id="one-past-end"),
    pytest.param(5000, 5000, "start=5000, end=5000", id="zero-width-past-end"),
    pytest.param(99999999, 5, "start=99999999", id="start-past-end-with-start-gt-end"),
    pytest.param(-1, 1, "start=-1", id="negative-start"),
    pytest.param(0, -1, "end=-1", id="negative-end"),
    pytest.param(
        4400, models.Range([4500, 4600]), "end=[4500, 4600]", id="definite-end-past-end"
    ),
    pytest.param(
        4400,
        models.Range([4561, None]),
        "end=[4561, None]",
        id="indefinite-end-lower-bound-past-end",
    ),
    pytest.param(
        models.Range([-5, 10]), 20, "start=[-5, 10]", id="range-with-negative-member"
    ),
]


@pytest.mark.parametrize(("start", "end"), LOCATION_BOUNDS_VALID)
def test_validate_location_bounds_valid(
    start: int | models.Range | None,
    end: int | models.Range | None,
) -> None:
    dp = _StubDataProxy({BOUNDS_SEQ_ID: BOUNDS_SEQ_LEN})
    dp.validate_location_bounds(BOUNDS_SEQ_ID, start, end)


@pytest.mark.parametrize(("start", "end", "detail"), LOCATION_BOUNDS_INVALID)
def test_validate_location_bounds_invalid(
    start: int | models.Range | None,
    end: int | models.Range | None,
    detail: str,
) -> None:
    dp = _StubDataProxy({BOUNDS_SEQ_ID: BOUNDS_SEQ_LEN})
    expected_msg = (
        f"Location out of bounds on {BOUNDS_SEQ_ID}: {detail} "
        f"not within [0, {BOUNDS_SEQ_LEN}]"
    )
    with pytest.raises(DataProxyValidationError, match=f"^{re.escape(expected_msg)}$"):
        dp.validate_location_bounds(BOUNDS_SEQ_ID, start, end)
