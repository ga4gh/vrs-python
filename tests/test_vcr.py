import pytest

vcr_test = (
    ("rest_dataproxy", "NC_000013.11"),
    ("dataproxy", "refseq:NC_000013.11")
)


@pytest.mark.parametrize("dp,identifier", vcr_test)
@pytest.mark.vcr
def test_vcrtest(dp, identifier, request):
    dp = request.getfixturevalue(dp)
    seq = dp.get_sequence(identifier, 50_000_000, 50_000_050)
    assert len(seq) == 50
    assert seq == "TTAGGTGTTTAGATGATTTCTAAGATGCTTTTAAGCCCAGTATTTCTATT"
