"""Check that forbidden operations fail at the call, including HGVS fallbacks."""

import socket

import psycopg2
import pytest
import requests

from ga4gh.vrs import models
from ga4gh.vrs.utils.hgvs_tools import HgvsTools


@pytest.mark.parametrize("client", ["http", "connect_ex", "postgres"])
def test_live_connections_fail(network_blocked: bool, client: str) -> None:
    if not network_blocked:
        pytest.skip("live connections are enabled for recording")
    message = "Live PostgreSQL" if client == "postgres" else "Live network"

    def connect_with_fallback() -> None:
        # An ordinary application error handler must not hide the failure.
        try:
            if client == "http":
                requests.get("http://127.0.0.1:1/", timeout=1)
            elif client == "connect_ex":
                with socket.socket() as sock:
                    sock.connect_ex(("127.0.0.1", 1))
            else:
                psycopg2.connect("postgresql://anonymous@127.0.0.1:1/uta")
        except Exception:
            return

    with pytest.raises(pytest.fail.Exception, match=message):
        connect_with_fallback()


def test_hgvs_cache_miss_fails_through_fallback(hgvs_cache_mode: str) -> None:
    if hgvs_cache_mode != "run":
        pytest.skip("HGVS recording uses live data")
    # An insertion needs no data-proxy fetch before reaching HGVS normalization.
    # This accession is deliberately absent from the recorded provider data.
    allele = models.Allele(
        location=models.SequenceLocation(start=10, end=10),
        state=models.LiteralSequenceExpression(sequence="A"),
    )
    tools = HgvsTools()
    with pytest.raises(pytest.fail.Exception, match="Missing HGVS recording"):
        tools._to_sequence_variant(allele, "c", "unused", "NM_000000.0")
