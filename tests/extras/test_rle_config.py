import pytest

from ga4gh.vrs import config


@pytest.mark.parametrize("limit", ["-1", "invalid", "1.5"])
def test_invalid_rle_sequence_limit_environment(limit, monkeypatch):
    monkeypatch.setenv("GA4GH_VRS_RLE_SEQ_LIMIT", limit)
    with pytest.raises(ValueError, match="GA4GH_VRS_RLE_SEQ_LIMIT"):
        config._get_rle_seq_limit()
