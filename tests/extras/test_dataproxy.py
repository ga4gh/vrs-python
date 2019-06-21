from vcr_support import vcr


@vcr.use_cassette()
def test_dataproxy_rest(dataproxy):
    r = dataproxy.get_metadata("NM_000551.3")
    assert 4560 == r["length"]
    assert   "VMC:GS_v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_" in r["aliases"]
    assert "ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_" in r["aliases"]

    r = dataproxy.get_metadata("NC_000013.11")
    assert 114364328 == r["length"]
    assert   "VMC:GS__0wi-qoDrvram155UmcSC-zA5ZK4fpLT" in r["aliases"]
    assert "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT" in r["aliases"]
    
