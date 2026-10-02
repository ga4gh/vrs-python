from importlib import import_module, reload
from pathlib import Path

import pysam
import pytest
import vcr

from ga4gh.vrs import config
from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy


@pytest.fixture
def configured_modules(monkeypatch, request):
    """Reload import-time defaults, restoring them after each environment case."""
    modules = [
        config,
        import_module("ga4gh.vrs.normalize"),
        import_module("ga4gh.vrs.extras.translator"),
        import_module("ga4gh.vrs.extras.annotator.vcf"),
    ]
    original_namespaces = [vars(module).copy() for module in modules]
    try:
        with monkeypatch.context() as env:
            env.delenv("GA4GH_VRS_RLE_SEQ_LIMIT", raising=False)
            if request.param is not None:
                env.setenv("GA4GH_VRS_RLE_SEQ_LIMIT", request.param)
            for module in modules:
                reload(module)
            yield modules[1:]
    finally:
        # Preserve class identities already imported by other test modules.
        for module, namespace in zip(modules, original_namespaces, strict=True):
            vars(module).clear()
            vars(module).update(namespace)


@pytest.mark.parametrize(
    ("configured_modules", "expected_sequence"),
    [
        (None, "CTTTCTTT"),
        ("0", None),
        ("4", None),
        ("8", "CTTTCTTT"),
        ("none", "CTTTCTTT"),
        ("NoNe", "CTTTCTTT"),
    ],
    indirect=["configured_modules"],
)
def test_rle_sequence_limit_environment(
    configured_modules, expected_sequence, tmp_path
):
    """Normalization, translation and VCF output share the configured default."""
    normalization, translator, vcf_module = configured_modules
    output = tmp_path / "annotated.vcf"
    with vcr.use_cassette(
        "tests/extras/cassettes/test_annotate_vcf_rle.yaml",
        record_mode="none",
        allow_playback_repeats=True,
    ):
        proxy = SeqRepoRESTDataProxy(
            base_url="http://localhost:5000/seqrepo",
            disable_healthcheck=True,
        )
        tlr = translator.AlleleTranslator(proxy)
        variant = "1-102995989-CTTT-CTTTCTTT"
        state = tlr.translate_from(variant, fmt="gnomad").state
        assert (
            state.sequence.root if state.sequence is not None else None
        ) == expected_sequence
        literal = tlr.translate_from(variant, fmt="gnomad", do_normalize=False)
        normalized = normalization.normalize(literal, proxy).state
        assert (
            normalized.sequence.root if normalized.sequence is not None else None
        ) == expected_sequence
        explicit = tlr.translate_from(variant, fmt="gnomad", rle_seq_limit=None).state
        assert explicit.sequence.root == "CTTTCTTT"
        vcf_module.VcfAnnotator(proxy).annotate(
            Path("tests/extras/data/test_rle.vcf"),
            output,
            vrs_attributes=True,
        )
        with pysam.VariantFile(output) as vcf:
            records = list(vcf)
            assert records[1].info["VRS_States"][-1] == (expected_sequence or ".")


@pytest.mark.parametrize("limit", ["-1", "invalid", "1.5"])
def test_invalid_rle_sequence_limit_environment(limit, monkeypatch):
    monkeypatch.setenv("GA4GH_VRS_RLE_SEQ_LIMIT", limit)
    with pytest.raises(ValueError, match="GA4GH_VRS_RLE_SEQ_LIMIT"):
        config._get_rle_seq_limit()
