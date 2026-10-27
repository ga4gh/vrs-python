"""Tests for RefgetSequenceDataProxy, CompositeDataProxy, and _extract_digest."""

import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock
from urllib.parse import urlparse

import pytest

from ga4gh.vrs.dataproxy import (
    CompositeDataProxy,
    DataProxyValidationError,
    create_dataproxy,
)
from ga4gh.vrs.dataproxy_refget import (
    RefgetAliasDataProxy,
    RefgetSequenceDataProxy,
    _create_refget_dataproxy,
    _extract_digest,
    _resolve_refget_cache_base_dir,
    _resolve_remote_cache_dir,
    _slug_remote_url,
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
        store.get_sequence_metadata.return_value = metadata
        store.get_sequence_metadata_by_alias.return_value = metadata
        store.get_substring.return_value = SEQ
    else:
        store.get_sequence.return_value = None
        store.get_sequence_metadata.return_value = None
        store.get_sequence_metadata_by_alias.return_value = None
        store.get_substring.return_value = None
    return store


def _install_fake_gtars(
    monkeypatch: pytest.MonkeyPatch,
    *,
    local_store: MagicMock | None = None,
    remote_store: MagicMock | None = None,
) -> tuple[MagicMock, MagicMock, type]:
    """Install a fake gtars.refget module for factory tests."""
    local_store = local_store or MagicMock()
    remote_store = remote_store or MagicMock()

    class FakeRefgetStore:
        open_local = MagicMock(return_value=local_store)
        open_remote = MagicMock(return_value=remote_store)

    fake_gtars = SimpleNamespace()
    fake_refget = SimpleNamespace(RefgetStore=FakeRefgetStore)
    monkeypatch.setitem(sys.modules, "gtars", fake_gtars)
    monkeypatch.setitem(sys.modules, "gtars.refget", fake_refget)
    return local_store, remote_store, FakeRefgetStore


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

    def test_empty_interval_returns_empty(self) -> None:
        """A zero-width interval, before or after clamping, is the empty string."""
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        assert dp.get_sequence(GA4GH_ID, 10, 10) == ""
        assert dp.get_sequence(GA4GH_ID, -5, 0) == ""
        store.get_substring.assert_not_called()

    def test_start_after_end_rejected(self) -> None:
        """A start after the end is rejected before the store is consulted."""
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        with pytest.raises(DataProxyValidationError, match="start=10, end=5"):
            dp.get_sequence(GA4GH_ID, 10, 5)
        store.get_substring.assert_not_called()
        store.get_sequence_metadata.assert_not_called()


# --- RefgetSequenceDataProxy.get_sequence_length ---


class TestRefgetGetSequenceLength:
    def test_found(self) -> None:
        store = _make_mock_store()
        dp = RefgetSequenceDataProxy(store)
        assert dp.get_sequence_length(GA4GH_ID) == SEQ_LENGTH
        store.get_sequence_metadata.assert_called_once_with(DIGEST)
        store.get_sequence.assert_not_called()

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


# --- RefgetAliasDataProxy ---


class TestRefgetAliasDataProxy:
    def test_resolve_alias_uses_metadata_only(self) -> None:
        store = _make_mock_store()
        dp = RefgetAliasDataProxy(store)

        assert dp._resolve_to_digest("refseq:NM_000551.3") == DIGEST
        store.get_sequence_metadata_by_alias.assert_called_once_with(
            "refseq", "NM_000551.3"
        )
        store.get_sequence_by_alias.assert_not_called()

    def test_derive_refget_accession_skips_reverse_alias_scan(self) -> None:
        store = _make_mock_store()
        dp = RefgetAliasDataProxy(store)

        assert dp.derive_refget_accession("refseq:NM_000551.3") == SQ_ID
        store.get_aliases_for_sequence.assert_not_called()

    @pytest.mark.parametrize(
        ("bare", "namespace"),
        [("NM_000551.3", "refseq"), ("LRG_1", "lrg"), ("LRG_1t1", "lrg")],
    )
    def test_bare_accession_namespace_is_inferred(
        self, bare: str, namespace: str
    ) -> None:
        store = _make_mock_store()
        dp = RefgetAliasDataProxy(store)

        assert dp.derive_refget_accession(bare) == SQ_ID
        store.get_sequence_metadata_by_alias.assert_called_once_with(namespace, bare)

    def test_bare_accession_with_no_inferable_namespace_is_a_miss(self) -> None:
        store = _make_mock_store()
        dp = RefgetAliasDataProxy(store)

        assert dp.derive_refget_accession("NOT_AN_ACCESSION") is None
        store.get_sequence_metadata_by_alias.assert_not_called()


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
    alias_proxy.derive_refget_accession.return_value = SQ_ID

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
        alias_proxy.derive_refget_accession.assert_called_once_with("NM_000551.3")

    def test_alias_not_found(self) -> None:
        """Non-GA4GH identifier with no alias raises KeyError."""
        composite, _, alias_proxy = _make_composite()
        alias_proxy.derive_refget_accession.return_value = None
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
        alias_proxy.derive_refget_accession.assert_called_once_with("NM_000551.3")


class TestCompositeGetAliases:
    def test_delegates_to_alias_proxy(self) -> None:
        composite, _, alias_proxy = _make_composite()
        aliases = composite.get_aliases("NM_000551.3")
        assert GA4GH_ID in aliases
        alias_proxy.get_aliases.assert_called_once_with("NM_000551.3")

    def test_derive_refget_accession_delegates_to_alias_proxy(self) -> None:
        composite, _, alias_proxy = _make_composite()

        assert composite.derive_refget_accession("NM_000551.3") == SQ_ID
        alias_proxy.derive_refget_accession.assert_called_once_with("NM_000551.3")
        alias_proxy.get_aliases.assert_not_called()


class TestRefgetRemoteCacheHelpers:
    def test_slug_remote_url_is_stable(self) -> None:
        url = "https://example.org/refget-store"
        assert _slug_remote_url(url) == _slug_remote_url(url)

    def test_slug_remote_url_distinguishes_urls(self) -> None:
        assert _slug_remote_url("https://example.org/store-a") != _slug_remote_url(
            "https://example.org/store-b"
        )

    def test_resolve_cache_base_uses_env_override(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        cache_root = tmp_path / "custom-cache"
        monkeypatch.setenv("GA4GH_VRS_REFGET_CACHE_DIR", str(cache_root))
        assert _resolve_refget_cache_base_dir() == cache_root

    def test_resolve_cache_base_uses_xdg(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        monkeypatch.delenv("GA4GH_VRS_REFGET_CACHE_DIR", raising=False)
        monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
        assert _resolve_refget_cache_base_dir() == tmp_path / "xdg" / "gtars" / "refget"

    def test_resolve_cache_base_falls_back_to_home(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        monkeypatch.delenv("GA4GH_VRS_REFGET_CACHE_DIR", raising=False)
        monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
        monkeypatch.setenv("HOME", str(tmp_path / "home"))
        assert (
            _resolve_refget_cache_base_dir()
            == tmp_path / "home" / ".cache" / "gtars" / "refget"
        )

    def test_resolve_remote_cache_dir_uses_base_plus_slug(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        cache_root = tmp_path / "cache-root"
        remote_url = "https://example.org/refget-store"
        monkeypatch.setenv("GA4GH_VRS_REFGET_CACHE_DIR", str(cache_root))
        assert _resolve_remote_cache_dir(remote_url) == cache_root / _slug_remote_url(
            remote_url
        )


class TestCreateRefgetDataProxyFactory:
    def test_file_uri_uses_open_local(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        local_store, _, fake_refget_store = _install_fake_gtars(monkeypatch)
        parsed_uri = urlparse(f"refget+file://{tmp_path}/store")

        proxy = _create_refget_dataproxy(parsed_uri, "file")

        assert isinstance(proxy, CompositeDataProxy)
        fake_refget_store.open_local.assert_called_once_with(parsed_uri.path)
        fake_refget_store.open_remote.assert_not_called()
        local_store.set_quiet.assert_called_once_with(True)
        local_store.pull_aliases.assert_called_once_with()

    def test_http_uri_uses_open_remote_with_env_cache_root(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        _, remote_store, fake_refget_store = _install_fake_gtars(monkeypatch)
        cache_root = tmp_path / "cache-root"
        remote_url = "http://example.org:8080/refget-store"
        parsed_uri = urlparse(f"refget+{remote_url}")
        monkeypatch.setenv("GA4GH_VRS_REFGET_CACHE_DIR", str(cache_root))

        proxy = _create_refget_dataproxy(parsed_uri, "http")

        assert isinstance(proxy, CompositeDataProxy)
        expected_cache_dir = cache_root / _slug_remote_url(remote_url)
        fake_refget_store.open_local.assert_not_called()
        fake_refget_store.open_remote.assert_called_once_with(
            str(expected_cache_dir), remote_url
        )
        assert expected_cache_dir.is_dir()
        remote_store.pull_aliases.assert_called_once_with()

    def test_https_uri_dispatches_via_create_dataproxy(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        _, remote_store, fake_refget_store = _install_fake_gtars(monkeypatch)
        cache_root = tmp_path / "cache-root"
        remote_url = "https://example.org/refget-store"
        monkeypatch.setenv("GA4GH_VRS_REFGET_CACHE_DIR", str(cache_root))

        proxy = create_dataproxy(f"refget+{remote_url}")

        assert isinstance(proxy, CompositeDataProxy)
        fake_refget_store.open_remote.assert_called_once_with(
            str(cache_root / _slug_remote_url(remote_url)),
            remote_url,
        )
        remote_store.pull_aliases.assert_called_once_with()

    def test_local_refget_uri_ignores_cache_env_var(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        local_store, _, fake_refget_store = _install_fake_gtars(monkeypatch)
        monkeypatch.setenv("GA4GH_VRS_REFGET_CACHE_DIR", str(tmp_path / "cache-root"))

        proxy = create_dataproxy(f"refget+file://{tmp_path}/store")

        assert isinstance(proxy, CompositeDataProxy)
        fake_refget_store.open_local.assert_called_once()
        fake_refget_store.open_remote.assert_not_called()
        local_store.pull_aliases.assert_called_once_with()
