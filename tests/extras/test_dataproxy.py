import pytest


@pytest.mark.parametrize("dp", ("rest_dataproxy","dataproxy"))
@pytest.mark.vcr
def test_data_proxies(dp, request):
    dataproxy = request.getfixturevalue(dp)
    r = dataproxy.get_metadata("NM_000551.3")
    assert 4560 == r["length"]
    assert "ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_" in r["aliases"]

    r = dataproxy.get_metadata("NC_000013.11")
    assert 114364328 == r["length"]
    assert "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT" in r["aliases"]

    seq = dataproxy.get_sequence("ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_")
    assert seq.startswith("CCTCGCCTCCGTTACAACGGCCTACGGTGCTGGAGGATCCTTCTGCGCAC")

    seq = dataproxy.get_sequence("ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_", start=0, end=50)
    assert seq == "CCTCGCCTCCGTTACAACGGCCTACGGTGCTGGAGGATCCTTCTGCGCAC"
