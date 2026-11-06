"""Tests for ga4gh.vrs.utils.hgvs_tools."""

from unittest.mock import Mock

import hgvs.parser
import pytest

from ga4gh.vrs.utils.hgvs_tools import HgvsTools


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


@pytest.fixture
def lifecycle(monkeypatch):
    """Mock the UTA connection, normalizer, and variant mapper."""
    connection = Mock()
    monkeypatch.setattr("hgvs.dataproviders.uta.connect", Mock(return_value=connection))
    monkeypatch.setattr("hgvs.normalizer.Normalizer", Mock())
    monkeypatch.setattr("hgvs.variantmapper.VariantMapper", Mock())
    return connection


def test_close_preserves_borrowed_resources(lifecycle):
    """Close UTA once, preserve references, and leave the borrowed proxy open."""
    proxy = Mock()
    tools = HgvsTools(proxy)
    normalizer, mapper = tools.normalizer, tools.variant_mapper
    tools.close()
    tools.close()
    lifecycle.close.assert_called_once_with()
    proxy.close.assert_not_called()
    assert tools.data_proxy is proxy
    assert tools.normalizer is normalizer
    assert tools.variant_mapper is mapper
    with pytest.raises(RuntimeError, match="closed"):
        tools.__enter__()


@pytest.mark.parametrize("fail", [False, True])
def test_context_cleanup(lifecycle, fail):
    """Context exit closes the owned UTA connection and propagates body errors."""
    tools = HgvsTools()
    from contextlib import nullcontext

    error = ValueError("body error")
    expectation = (
        pytest.raises(ValueError, match="body error") if fail else nullcontext()
    )
    with expectation, tools as entered:
        assert entered is tools
        if fail:
            raise error
    lifecycle.close.assert_called_once_with()


@pytest.mark.parametrize("collaborator", ["Normalizer", "VariantMapper"])
def test_initialization_failure_closes_connection(lifecycle, monkeypatch, collaborator):
    """Close the owned UTA connection if normalizer or variant mapper setup fails."""
    module = "normalizer" if collaborator == "Normalizer" else "variantmapper"
    monkeypatch.setattr(
        f"hgvs.{module}.{collaborator}", Mock(side_effect=ValueError("init"))
    )
    with pytest.raises(ValueError, match="init"):
        HgvsTools()
    lifecycle.close.assert_called_once_with()


def test_failed_cleanup_can_retry(lifecycle):
    """Propagate cleanup errors and allow another close attempt."""
    tools = HgvsTools()
    lifecycle.close.side_effect = [RuntimeError("cleanup"), None]
    with pytest.raises(RuntimeError, match="cleanup"):
        tools.close()
    tools.close()
    tools.close()
    assert lifecycle.close.call_count == 2
