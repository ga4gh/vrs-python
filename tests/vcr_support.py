import os
from pathlib import Path

import vcr as vcrpy

test_dir = Path(__file__).parent
test_data_dir = Path(test_dir) / "data" / "cassettes"

vcr = vcrpy.VCR(
    cassette_library_dir=test_data_dir,
    record_mode=os.environ.get("VCR_RECORD_MODE", "new_episodes"),
    decode_compressed_response=True,
)
