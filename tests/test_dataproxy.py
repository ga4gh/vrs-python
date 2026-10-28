import re
import sys
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import (
    AliasDataProxy,
    DataProxy,
    DataProxyValidationError,
    UnresolvableAccessionError,
    coerce_accession_namespace,
    create_dataproxy,
    require_refget_accession,
    resolve_cache_size,
)
from ga4gh.vrs.dataproxy_seqrepo import SeqRepoDataProxy, SeqRepoRESTDataProxy


@pytest.mark.parametrize(
    ("ac", "expected"),
    [
        ("LRG_1", "lrg:LRG_1"),
        ("LRG_1g", "lrg:LRG_1g"),
        ("LRG_199t1", "lrg:LRG_199t1"),
        ("LRG_199p1", "lrg:LRG_199p1"),
        ("NM_000551.3", "refseq:NM_000551.3"),
        ("ENST00000530893.6", "ensembl:ENST00000530893.6"),
        ("lrg:LRG_1t1", "lrg:LRG_1t1"),
    ],
)
def test_coerce_accession_namespace(ac: str, expected: str) -> None:
    assert coerce_accession_namespace(ac) == expected


@pytest.mark.parametrize(
    "ac", ["LRG_", "LRG_1x1", "LRG_1t", "LRG1", "NOT_AN_ACCESSION"]
)
def test_coerce_accession_namespace_rejects(ac: str) -> None:
    with pytest.raises(ValueError, match="Could not infer namespace"):
        coerce_accession_namespace(ac)


@pytest.mark.parametrize("dp", ["rest_dataproxy", "dataproxy"])
@pytest.mark.vcr
def test_data_proxies(dp, request):
    dataproxy = request.getfixturevalue(dp)

    length = dataproxy.get_sequence_length("NM_000551.3")
    assert length == 4560

    aliases = dataproxy.get_aliases("NM_000551.3")
    assert "ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_" in aliases

    length = dataproxy.get_sequence_length("NC_000013.11")
    assert length == 114364328

    aliases = dataproxy.get_aliases("NC_000013.11")
    assert "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT" in aliases

    seq = dataproxy.get_sequence("ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_")
    assert seq.startswith("CCTCGCCTCCGTTACAACGGCCTACGGTGCTGGAGGATCCTTCTGCGCAC")

    seq = dataproxy.get_sequence(
        "ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_", start=0, end=50
    )
    assert seq == "CCTCGCCTCCGTTACAACGGCCTACGGTGCTGGAGGATCCTTCTGCGCAC"


def test_data_proxy_configs(monkeypatch):
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

    monkeypatch.delenv("GA4GH_VRS_DATAPROXY_URI", raising=False)
    with pytest.raises(
        ValueError,
        match="No data proxy URI provided or found in GA4GH_VRS_DATAPROXY_URI",
    ):
        create_dataproxy(None)

    # check that fallback on env var works
    monkeypatch.setenv("GA4GH_VRS_DATAPROXY_URI", "seqrepo+:tests/data/seqrepo/latest")
    create_dataproxy(None)

    # check that method arg takes precedence by passing an invalid arg
    with pytest.raises(
        ValueError,
        match=re.escape(
            "create_dataproxy scheme must include provider (e.g., `seqrepo+http:...`)"
        ),
    ):
        create_dataproxy("file:///path/to/seqrepo/root")


def test_refget_remote_uri_configs(monkeypatch, tmp_path):
    remote_store = MagicMock()

    class FakeRefgetStore:
        open_local = MagicMock()
        open_remote = MagicMock(return_value=remote_store)

    monkeypatch.setitem(sys.modules, "gtars", SimpleNamespace())
    monkeypatch.setitem(
        sys.modules,
        "gtars.refget",
        SimpleNamespace(RefgetStore=FakeRefgetStore),
    )
    monkeypatch.setenv("GA4GH_VRS_REFGET_CACHE_DIR", str(tmp_path / "cache-root"))

    create_dataproxy("refget+http://example.org/refget-store")
    create_dataproxy("refget+https://example.org/refget-store")

    assert FakeRefgetStore.open_remote.call_count == 2

    with pytest.raises(
        ValueError,
        match=re.escape("Refget URI protocol 'fake-scheme' not supported"),
    ):
        create_dataproxy("refget+fake-scheme://example.org/refget-store")


class TestResolveCacheSize:
    """Cache sizes are resolved from the environment at import time."""

    ENV_VAR = "GA4GH_VRS_TEST_CACHE_SIZE"

    @pytest.mark.parametrize(
        ("env_value", "expected"),
        [
            (None, 8192),  # unset falls back to the default
            ("", 8192),  # empty is treated as unset
            ("512", 512),
            ("0", 0),  # functools semantics: caching disabled
            ("-1", None),  # unbounded
            ("-99", None),
            ("none", None),
            ("NONE", None),
            ("  none  ", None),
            ("  512  ", 512),
        ],
    )
    def test_valid_values(
        self,
        monkeypatch: pytest.MonkeyPatch,
        env_value: str | None,
        expected: int | None,
    ) -> None:
        if env_value is None:
            monkeypatch.delenv(self.ENV_VAR, raising=False)
        else:
            monkeypatch.setenv(self.ENV_VAR, env_value)
        assert resolve_cache_size(self.ENV_VAR, "8192") == expected

    @pytest.mark.parametrize("env_value", ["abc", "1.5", "8k"])
    def test_invalid_values_raise(
        self, monkeypatch: pytest.MonkeyPatch, env_value: str
    ) -> None:
        monkeypatch.setenv(self.ENV_VAR, env_value)
        with pytest.raises(ValueError, match="must be an integer or 'none'"):
            resolve_cache_size(self.ENV_VAR, "8192")


class TestDeriveRefgetAccession:
    """A sequence has exactly one refget accession, so >1 ga4gh alias is a bug."""

    SQ_A = "ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_"
    SQ_B = "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"

    @staticmethod
    def _proxy(aliases: list[str]) -> AliasDataProxy:
        class _Proxy(AliasDataProxy):
            def get_aliases(self, identifier: str) -> list[str]:  # noqa: ARG002
                return aliases

        return _Proxy()

    @pytest.mark.parametrize(
        ("aliases", "expected"),
        [
            ([SQ_A, "refseq:NM_000551.3"], "SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_"),
            (["refseq:NM_000551.3"], None),  # no ga4gh alias is a normal miss
        ],
    )
    def test_resolves(self, aliases: list[str], expected: str | None) -> None:
        proxy = self._proxy(aliases)
        assert proxy.derive_refget_accession("refseq:NM_000551.3") == expected

    def test_multiple_ga4gh_aliases_raises(self) -> None:
        proxy = self._proxy([self.SQ_A, self.SQ_B, "refseq:NM_000551.3"])
        with pytest.raises(ValueError, match="Expected at most one ga4gh alias"):
            proxy.derive_refget_accession("refseq:NM_000551.3")

    def test_uninferrable_namespace_is_a_miss(self) -> None:
        proxy = self._proxy([self.SQ_A])
        assert proxy.derive_refget_accession("NOT_AN_ACCESSION") is None

    def test_require_refget_accession(self) -> None:
        assert (
            require_refget_accession(self._proxy([self.SQ_A]), "refseq:NM_000551.3")
            == "SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_"
        )
        with pytest.raises(UnresolvableAccessionError) as excinfo:
            require_refget_accession(self._proxy([]), "refseq:NM_000551.3")
        assert isinstance(excinfo.value, KeyError)
        assert excinfo.value.accession == "refseq:NM_000551.3"
        assert str(excinfo.value) == (
            "Unable to resolve refseq:NM_000551.3 to a refget accession"
        )


class _StubDataProxy(DataProxy):
    """Dataproxy serving only lengths and aliases, keyed by identifier"""

    def __init__(self, metadata: dict[str, dict]) -> None:
        self.metadata = metadata

    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        raise NotImplementedError

    def get_sequence_length(self, identifier: str) -> int:
        return self.metadata[identifier]["length"]

    def get_aliases(self, identifier: str) -> list[str]:
        return self.metadata[identifier]["aliases"]


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
    {
        "id": "overlapping-ranges",
        "start": models.Range([100, 200]),
        "end": models.Range([150, 300]),
    },
    {
        "id": "open-ranges-never-start-after-end",
        "start": models.Range([None, 4000]),
        "end": models.Range([10, None]),
    },
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


# In-bounds locations whose start is after their end. A Range is reported by its
# lowest possible start and highest possible end.
location_start_after_end_cases = [
    {"id": "start-after-end", "start": 4000, "end": 5, "detail": "start=4000, end=5"},
    {
        "id": "definite-ranges-disjoint",
        "start": models.Range([4000, 4100]),
        "end": models.Range([10, 20]),
        "detail": "start=4000, end=20",
    },
    {
        "id": "indefinite-ranges-disjoint",
        "start": models.Range([4000, None]),
        "end": models.Range([None, 20]),
        "detail": "start=4000, end=20",
    },
]


def _start_after_end_msg(identifier: str, detail: str) -> str:
    return f"Interval start must not be after end on {identifier}: {detail}"


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


@pytest.mark.parametrize("case", location_start_after_end_cases, ids=lambda c: c["id"])
def test_validate_location_bounds_start_after_end(case: dict) -> None:
    dp = _bounds_dp([BOUNDS_SEQ_ID, BOUNDS_REFGET_ID])
    expected_msg = _start_after_end_msg(BOUNDS_SEQ_ID, case["detail"])
    with pytest.raises(DataProxyValidationError, match=f"^{re.escape(expected_msg)}$"):
        dp.validate_location_bounds(BOUNDS_SEQ_ID, case["start"], case["end"])


def test_validate_ref_seq_checks_bounds_before_fetching() -> None:
    # The stub's get_sequence raises, so this also shows nothing is fetched
    dp = _bounds_dp([BOUNDS_SEQ_ID, BOUNDS_REFGET_ID])
    expected_msg = (
        f"Location out of bounds on {BOUNDS_SEQ_ID} ({BOUNDS_REFGET_ID}): "
        f"end=4561 not within [0, {BOUNDS_SEQ_LEN}]"
    )
    with pytest.raises(DataProxyValidationError, match=f"^{re.escape(expected_msg)}$"):
        dp.validate_ref_seq(BOUNDS_SEQ_ID, 4560, 4561, "A", require_validation=False)


def test_seqrepo_get_sequence_rejects_start_after_end_before_the_backend() -> None:
    sr = MagicMock()
    expected_msg = _start_after_end_msg(BOUNDS_SEQ_ID, "start=4000, end=5")
    with pytest.raises(DataProxyValidationError, match=f"^{re.escape(expected_msg)}$"):
        SeqRepoDataProxy(sr).get_sequence(BOUNDS_SEQ_ID, 4000, 5)
    sr.fetch_uri.assert_not_called()


def test_seqrepo_rest_get_sequence_rejects_start_after_end_before_the_backend(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    get = MagicMock()
    monkeypatch.setattr("ga4gh.vrs.dataproxy_seqrepo.requests.get", get)
    dp = SeqRepoRESTDataProxy("http://localhost:5000/seqrepo", disable_healthcheck=True)
    expected_msg = _start_after_end_msg(BOUNDS_SEQ_ID, "start=4000, end=5")
    with pytest.raises(DataProxyValidationError, match=f"^{re.escape(expected_msg)}$"):
        dp.get_sequence(BOUNDS_SEQ_ID, 4000, 5)
    get.assert_not_called()


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
