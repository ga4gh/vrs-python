"""Process-level configuration for VRS-Python."""

import os


def _get_rle_seq_limit() -> int | None:
    """Read the maximum optional RLE sequence length from the environment."""
    value = os.environ.get("GA4GH_VRS_RLE_SEQ_LIMIT", "50")
    if value.lower() == "none":
        return None
    message = "GA4GH_VRS_RLE_SEQ_LIMIT must be a non-negative integer or 'none'"
    try:
        limit = int(value)
    except ValueError as exc:
        raise ValueError(message) from exc
    if limit < 0:
        raise ValueError(message)
    return limit


RLE_SEQ_LIMIT = _get_rle_seq_limit()
