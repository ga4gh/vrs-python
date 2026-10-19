"""Tests for ga4gh.vrs.utils.hgvs_tools."""

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
