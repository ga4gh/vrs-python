"""Tests for ga4gh.vrs.utils.hgvs_tools."""

import os

import hgvs.dataproviders.uta
import hgvs.parser
import pytest

from ga4gh.vrs.utils.hgvs_tools import HgvsTools, _connect_with_cache


@pytest.fixture(scope="module")
def hgvs_tools():
    # is_intronic only parses expressions, so build HgvsTools without the
    # network-dependent __init__ (which opens a UTA connection).
    tools = HgvsTools.__new__(HgvsTools)
    tools.parser = hgvs.parser.Parser()
    return tools


@pytest.mark.parametrize(
    ("hgvs_expr", "expected"),
    [
        # coding (c.)
        ("NM_000000.1:c.76-1G>A", True),  # intronic (acceptor side)
        ("NM_000000.1:c.87+1A>T", True),  # intronic (donor side)
        ("NM_000000.1:c.100A>T", False),  # exonic
        # RNA (r.) - previously missed: the interval is a plain Interval, not a
        # BaseOffsetInterval, so the old isinstance check returned False
        ("NR_000000.1:r.76-1g>a", True),  # intronic (acceptor side)
        ("NR_000000.1:r.77+1g>a", True),  # intronic (donor side)
        ("NR_000000.1:r.100a>u", False),  # exonic
        # RNA (r.) exonic - real-world ClinVar pins from issue #628 (offset 0)
        ("NR_001566.1:r.398_399del", False),
        ("NR_001566.1:r.245del", False),
        ("NM_001374385.1:r.2843_2931del", False),
        ("NM_001323289.2:r.2632c>a", False),
        # genomic (g.) - no base-offset positions, never intronic
        ("NC_000001.11:g.100A>T", False),
    ],
)
def test_is_intronic(hgvs_tools, hgvs_expr, expected):
    sv = hgvs_tools.parse(hgvs_expr)
    assert sv is not None
    assert hgvs_tools.is_intronic(sv) is expected


def test_hgvs_cache_run_mode_does_not_connect(monkeypatch):
    """With the tests' default hgvs cache (run mode), lookups are served from the cache
    file without connecting to UTA
    """
    if os.environ.get("VRS_HGVS_CACHE_MODE") != "run":
        pytest.skip("only applies when the hgvs cache is in run mode")

    def fail_connect(self):  # noqa: ARG001
        msg = "UTA connection attempted in hgvs cache run mode"
        raise AssertionError(msg)

    monkeypatch.setattr(hgvs.dataproviders.uta.UTA_postgresql, "_connect", fail_connect)
    _connect_with_cache.cache_clear()

    tools = HgvsTools()
    assert tools.uta_conn.get_tx_identity_info("NM_181798.1")["tx_ac"] == "NM_181798.1"


def test_hgvs_cache_data_provider_is_shared():
    """HgvsTools instances using the same cache share one data provider, so in learn
    mode they cannot overwrite each other's cache entries
    """
    assert HgvsTools().uta_conn is HgvsTools().uta_conn


hgvs_cache_env_invalid_cases = [
    {
        "id": "unknown-mode",
        "env": {"VRS_HGVS_CACHE_MODE": "lean"},
        "match": "must be one of",
    },
    {
        "id": "mode-without-file",
        "env": {"VRS_HGVS_CACHE_FILE": ""},
        "match": "VRS_HGVS_CACHE_FILE must be set",
    },
]


@pytest.mark.parametrize("case", hgvs_cache_env_invalid_cases, ids=lambda c: c["id"])
def test_hgvs_cache_env_invalid(monkeypatch, case):
    """Invalid hgvs cache settings raise rather than silently connecting to UTA"""
    for name, value in case["env"].items():
        monkeypatch.setenv(name, value)
    with pytest.raises(ValueError, match=case["match"]):
        HgvsTools()
