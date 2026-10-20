"""Tests for the session's offline guard (tests/offline_guard.py)."""

from typing import Any

import hgvs.dataproviders.uta
import psycopg2
import pytest
import requests
from hgvs.exceptions import HGVSDataNotAvailableError
from offline_guard import OfflineError, OfflineGuard

offline_guard_cases: list[dict[str, Any]] = [
    {
        "id": "http-request",
        "active_when": "blocking_network",
        "call": lambda _provider: requests.get("http://127.0.0.1:1/", timeout=1),
        "raises": OfflineError,
        "violation": "blocked: socket.connect(('127.0.0.1', 1))",
    },
    {
        "id": "postgres-connect",
        "active_when": "blocking_network",
        "call": lambda _provider: psycopg2.connect(
            "postgresql://anonymous@127.0.0.1:1/uta"
        ),
        "raises": OfflineError,
        "violation": "blocked: psycopg2.connect('postgresql://anonymous@127.0.0.1:1/uta')",
    },
    {
        "id": "hgvs-cache-miss",
        "active_when": "recording_hgvs_cache_misses",
        "call": lambda provider: provider.get_tx_info(
            "NM_000000.0", "NC_000001.11", "splign"
        ),
        "raises": HGVSDataNotAvailableError,
        "violation": "hgvs cache miss: get_tx_info('NM_000000.0', 'NC_000001.11', 'splign')",
    },
]


@pytest.mark.parametrize("case", offline_guard_cases, ids=lambda c: c["id"])
def test_offline_guard(
    case: dict[str, Any],
    offline_guard: OfflineGuard,
    hgvs_cached_data_provider: hgvs.dataproviders.uta.UTABase | None,
) -> None:
    """The offline guard blocks network connections and records hgvs cache misses, each
    as one violation
    """
    if not getattr(offline_guard, case["active_when"]):
        pytest.skip(f"the offline guard isn't {case['active_when'].replace('_', ' ')}")
    start = len(offline_guard.violations)
    with pytest.raises(case["raises"]):
        case["call"](hgvs_cached_data_provider)
    assert offline_guard.violations[start:] == [case["violation"]]
    # these violations are expected; keep them from failing this test at teardown
    del offline_guard.violations[start:]
