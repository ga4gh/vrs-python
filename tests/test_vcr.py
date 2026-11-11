import psycopg2
import pytest
import requests


@pytest.mark.vcr
def test_vcrtest(rest_dataproxy):
    seq = rest_dataproxy.get_sequence("NC_000013.11", 50_000_000, 50_000_050)
    assert len(seq) == 50
    assert seq == "TTAGGTGTTTAGATGATTTCTAAGATGCTTTTAAGCCCAGTATTTCTATT"


@pytest.mark.parametrize(
    "connect",
    [
        pytest.param(lambda: requests.get("http://127.0.0.1:1/", timeout=1), id="http"),
        pytest.param(
            lambda: psycopg2.connect(host="127.0.0.1", port=1, connect_timeout=1),
            id="postgres",
        ),
    ],
)
def test_live_connections_blocked(block_live_connections: bool, connect) -> None:
    if not block_live_connections:
        pytest.skip("live connections are enabled")
    with pytest.raises(pytest.fail.Exception):
        connect()
