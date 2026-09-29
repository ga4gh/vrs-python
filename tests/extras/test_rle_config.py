import json
import os
import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    ("limit", "expected_sequence"),
    [
        (None, "CTTTCTTT"),
        ("0", None),
        ("4", None),
        ("8", "CTTTCTTT"),
        ("none", "CTTTCTTT"),
    ],
)
def test_rle_sequence_limit_environment(limit, expected_sequence, tmp_path):
    """The translator and VCF writer honor the same process-level limit."""
    env = os.environ.copy()
    env.pop("GA4GH_VRS_RLE_SEQ_LIMIT", None)
    if limit is not None:
        env["GA4GH_VRS_RLE_SEQ_LIMIT"] = limit
    output = tmp_path / "annotated.vcf"
    result = subprocess.run(  # noqa: S603
        [
            sys.executable,
            "-c",
            """
import json
import sys
import pysam
import vcr
from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy
from ga4gh.vrs import normalize
from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.extras.annotator.vcf import VcfAnnotator
from pathlib import Path

with vcr.use_cassette(
    'tests/extras/cassettes/test_annotate_vcf_rle.yaml',
    record_mode='none', allow_playback_repeats=True,
):
    proxy = SeqRepoRESTDataProxy(base_url='http://localhost:5000/seqrepo', disable_healthcheck=True)
    variant = '1-102995989-CTTT-CTTTCTTT'
    state = AlleleTranslator(proxy).translate_from(variant, fmt='gnomad').state
    literal = AlleleTranslator(proxy).translate_from(variant, fmt='gnomad', do_normalize=False)
    normalized = normalize(literal, proxy).state
    explicit = AlleleTranslator(proxy).translate_from(variant, fmt='gnomad', rle_seq_limit=None).state
    VcfAnnotator(proxy).annotate(
        Path('tests/extras/data/test_rle.vcf'), Path(sys.argv[1]), vrs_attributes=True,
    )
    with pysam.VariantFile(sys.argv[1]) as vcf:
        records = list(vcf)
        print(json.dumps([
            state.sequence.root if state.sequence is not None else None,
            records[1].info['VRS_States'][-1],
            explicit.sequence.root,
            normalized.sequence.root if normalized.sequence is not None else None,
        ]))
""",
            str(output),
        ],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == [
        expected_sequence,
        expected_sequence or ".",
        "CTTTCTTT",
        expected_sequence,
    ]


@pytest.mark.parametrize("limit", ["-1", "invalid", "1.5"])
def test_invalid_rle_sequence_limit_environment(limit):
    env = dict(os.environ, GA4GH_VRS_RLE_SEQ_LIMIT=limit)
    result = subprocess.run(  # noqa: S603
        [sys.executable, "-c", "import ga4gh.vrs"],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "GA4GH_VRS_RLE_SEQ_LIMIT" in result.stderr
