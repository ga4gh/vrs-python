import re
import sys
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from ga4gh.vrs.dataproxy import (
    AliasDataProxy,
    create_dataproxy,
    resolve_cache_size,
)


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
