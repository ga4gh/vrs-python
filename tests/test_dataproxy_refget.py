"""Tests for RefgetSequenceDataProxy, CompositeDataProxy, and _extract_digest."""

from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from ga4gh.vrs.dataproxy import CompositeDataProxy
from ga4gh.vrs.dataproxy_refget import (
    RefgetSequenceDataProxy,
    _extract_digest,
)

# --- _extract_digest ---


@pytest.mark.parametrize(
    ("identifier", "expected"),
    [
        (
            "ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_",
            "v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_",
        ),
        ("SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_", "v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_"),
        ("ga4gh:SQ.", ""),
        ("SQ.", ""),
        ("NM_000551.3", None),
        ("refseq:NM_000551.3", None),
        ("MD5:215137b1973c1a5afcf86be7d999574a", None),
        ("", None),
        ("ga4gh:VA.something", None),
    ],
)
def test_extract_digest(identifier: str, expected: str | None) -> None:
    assert _extract_digest(identifier) == expected


# --- Mock helpers ---


DIGEST = "v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_"
GA4GH_ID = f"ga4gh:SQ.{DIGEST}"
SQ_ID = f"SQ.{DIGEST}"
SEQ = "CCTCGCCTCC"
SEQ_LENGTH = 4560


def _make_mock_store(
    *, has_sequence: bool = True, decode_returns: str | None = SEQ
) -> MagicMock:
    """Create a mock RefgetStore with configurable behavior."""
    store = MagicMock()
    if has_sequence:
        metadata = SimpleNamespace(length=SEQ_LENGTH, sha512t24u=DIGEST, md5="abc123")
        record = MagicMock()
        record.metadata = metadata
        record.decode.return_value = decode_returns
        store.get_sequence.return_value = record
        store.get_substring.return_value = SEQ
    else:
        store.get_sequence.return_value = None
        store.get_substring.return_value = None
    return store


# --- RefgetSequenceDataProxy.get_sequence ---


class TestRefgetGetSequence:
    def test_both_bounds(self) -> None:
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        result = dp.get_sequence(GA4GH_ID, 0, 10)
        assert result == SEQ
        store.get_substring.assert_called_once_with(DIGEST, 0, 10)

    def test_no_bounds(self) -> None:
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        result = dp.get_sequence(GA4GH_ID)
        assert result == SEQ
        store.get_sequence.assert_called_once_with(DIGEST)

    def test_start_only(self) -> None:
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        result = dp.get_sequence(GA4GH_ID, start=10)
        assert result == SEQ
        store.get_substring.assert_called_once_with(DIGEST, 10, SEQ_LENGTH)

    def test_end_only(self) -> None:
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        result = dp.get_sequence(GA4GH_ID, end=50)
        assert result == SEQ
        store.get_substring.assert_called_once_with(DIGEST, 0, 50)

    def test_sq_prefix(self) -> None:
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        result = dp.get_sequence(SQ_ID, 0, 10)
        assert result == SEQ
        store.get_substring.assert_called_once_with(DIGEST, 0, 10)

    def test_not_found_both_bounds(self) -> None:
        store = _make_mock_store(has_sequence=False)
        dp = RefgetSequenceDataProxy(store)
        with pytest.raises(KeyError, match=GA4GH_ID):
            dp.get_sequence(GA4GH_ID, 0, 10)

    def test_not_found_no_bounds(self) -> None:
        store = _make_mock_store(has_sequence=False)
        dp = RefgetSequenceDataProxy(store)
        with pytest.raises(KeyError, match=GA4GH_ID):
            dp.get_sequence(GA4GH_ID)

    def test_non_ga4gh_identifier(self) -> None:
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        with pytest.raises(KeyError, match="NM_000551.3"):
            dp.get_sequence("NM_000551.3")

    def test_decode_returns_none(self) -> None:
        store = _make_mock_store(decode_returns=None)
        dp = RefgetSequenceDataProxy(store)
        with pytest.raises(KeyError, match=GA4GH_ID):
            dp.get_sequence(GA4GH_ID)

    def test_negative_start_clamped(self) -> None:
        """Negative start is clamped to 0 (matching Python slice semantics)."""
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        result = dp.get_sequence(GA4GH_ID, -10, 5)
        assert result == SEQ
        store.get_substring.assert_called_once_with(DIGEST, 0, 5)

    def test_start_ge_end_returns_empty(self) -> None:
        """When clamped start >= end, return empty string."""
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        assert dp.get_sequence(GA4GH_ID, 10, 10) == ""
        assert dp.get_sequence(GA4GH_ID, -5, 0) == ""
        store.get_substring.assert_not_called()


# --- RefgetSequenceDataProxy.get_sequence_length ---


class TestRefgetGetSequenceLength:
    def test_found(self) -> None:
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        assert dp.get_sequence_length(GA4GH_ID) == SEQ_LENGTH

    def test_not_found(self) -> None:
        store = _make_mock_store(has_sequence=False)
        dp = RefgetSequenceDataProxy(store)
        with pytest.raises(KeyError, match=GA4GH_ID):
            dp.get_sequence_length(GA4GH_ID)

    def test_non_ga4gh_identifier(self) -> None:
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        with pytest.raises(KeyError, match="NM_000551.3"):
            dp.get_sequence_length("NM_000551.3")


# --- CompositeDataProxy ---


def _make_composite(
    *,
    seq_has_sequence: bool = True,
) -> tuple[CompositeDataProxy, MagicMock, MagicMock]:
    """Create a CompositeDataProxy with mock backends."""
    seq_store = _make_mock_store(has_sequence=seq_has_sequence)
    seq_proxy = RefgetSequenceDataProxy(seq_store)

    alias_proxy = MagicMock()
    alias_proxy.get_aliases.return_value = [
        "refseq:NM_000551.3",
        GA4GH_ID,
        "MD5:abc123",
    ]
    alias_proxy.translate_sequence_identifier.return_value = [GA4GH_ID]

    composite = CompositeDataProxy(seq_proxy, alias_proxy)
    return composite, seq_store, alias_proxy


class TestCompositeGetSequence:
    def test_direct_ga4gh_hit(self) -> None:
        """GA4GH identifier succeeds on first try, no alias lookup."""
        composite, seq_store, alias_proxy = _make_composite()
        result = composite.get_sequence(GA4GH_ID, 0, 10)
        assert result == SEQ
        alias_proxy.translate_sequence_identifier.assert_not_called()

    def test_alias_fallback(self) -> None:
        """Non-GA4GH identifier triggers alias resolution then succeeds."""
        composite, seq_store, alias_proxy = _make_composite()
        result = composite.get_sequence("NM_000551.3", 0, 10)
        assert result == SEQ
        alias_proxy.translate_sequence_identifier.assert_called_once_with(
            "NM_000551.3", namespace="ga4gh"
        )

    def test_alias_not_found(self) -> None:
        """Non-GA4GH identifier with no alias raises KeyError."""
        composite, _, alias_proxy = _make_composite()
        alias_proxy.translate_sequence_identifier.return_value = []
        with pytest.raises(KeyError, match="NM_999999.1"):
            composite.get_sequence("NM_999999.1", 0, 10)


class TestCompositeGetSequenceLength:
    def test_direct_ga4gh_hit(self) -> None:
        composite, _, alias_proxy = _make_composite()
        assert composite.get_sequence_length(GA4GH_ID) == SEQ_LENGTH
        alias_proxy.translate_sequence_identifier.assert_not_called()

    def test_alias_fallback(self) -> None:
        composite, _, alias_proxy = _make_composite()
        assert composite.get_sequence_length("NM_000551.3") == SEQ_LENGTH
        alias_proxy.translate_sequence_identifier.assert_called_once_with(
            "NM_000551.3", namespace="ga4gh"
        )


class TestCompositeGetAliases:
    def test_delegates_to_alias_proxy(self) -> None:
        composite, _, alias_proxy = _make_composite()
        aliases = composite.get_aliases("NM_000551.3")
        assert GA4GH_ID in aliases
        alias_proxy.get_aliases.assert_called_once_with("NM_000551.3")
