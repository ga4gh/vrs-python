import os

import vcr as vcrpy

test_dir = os.path.dirname(__file__)
test_data_dir = os.path.join(test_dir, "data", "cassettes")

vcr = vcrpy.VCR(
    cassette_library_dir=test_data_dir,
    record_mode=os.environ.get("VCR_RECORD_MODE", "new_episodes"),
)
