def test_dataproxy_rest(dataproxy):
    r = dataproxy.get_metadata("NM_000551.3")
    assert 4560 == r["length"]
    assert 'VMC:GS_v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_' in r["aliases"]
    assert 'ga4gh:SQ/v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_' in r["aliases"]
