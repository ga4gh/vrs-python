"""Ensure proper functionality of VCFAnnotator"""

import gzip
import logging
import os
import re
from pathlib import Path

import pysam
import pytest

from ga4gh.vrs.dataproxy import DataProxyValidationError, SeqRepoRESTDataProxy
from ga4gh.vrs.extras.annotator.vcf import VcfAnnotator, VcfAnnotatorError

TEST_DATA_DIR = Path("tests/extras/data")


@pytest.fixture
def rest_dataproxy_fn_scope():
    """REST dataproxy scoped to individual test functions, rather than the entire session"""
    return SeqRepoRESTDataProxy(
        base_url=os.environ.get("SEQREPO_REST_URL", "http://localhost:5000/seqrepo")
    )


@pytest.fixture
def vcf_annotator(rest_dataproxy_fn_scope: SeqRepoRESTDataProxy):
    return VcfAnnotator(rest_dataproxy_fn_scope)


@pytest.fixture(scope="session")
def input_vcf():
    """Provide fixture for sample input VCF"""
    return TEST_DATA_DIR / "test_vcf_input.vcf"


def compare_vcfs(actual_vcf_path: Path, expected_vcf_path: Path):
    """Normalize version fields that change per build and compare the remaining content
    of two VCF files line-by-line.
    """
    version_section_pattern = re.compile(
        r"\[VRS version=[^;\]]+;VRS-Python version=[^\]]+\]"
    )

    def _mask_version_fields(line: str) -> str:
        """Replace the bracketed version section with a stable placeholder."""
        if not line.startswith("##INFO=<ID=VRS_Allele_IDs"):
            return line
        if not version_section_pattern.search(line):
            msg = f"Expected [VRS version=...;VRS-Python version=...] section in INFO header line {line}"
            raise AssertionError(msg)
        return version_section_pattern.sub(
            "[VRS version=<ignored>;VRS-Python version=<ignored>]", line
        )

    # Handle both gzipped and uncompressed VCF files
    def _open(path: Path):
        if path.suffix == ".gz":
            return gzip.open(path, "rt")
        return path.open("r")

    with _open(actual_vcf_path) as out_vcf, _open(expected_vcf_path) as expected_output:
        out_vcf_lines = out_vcf.readlines()
        expected_output_lines = expected_output.readlines()

    for actual_line, expected_line in zip(
        out_vcf_lines, expected_output_lines, strict=False
    ):
        actual_line = _mask_version_fields(actual_line)
        expected_line = _mask_version_fields(expected_line)
        assert actual_line == expected_line


@pytest.mark.vcr
def test_annotate_vcf_grch38_noattrs(
    vcf_annotator: VcfAnnotator, input_vcf: Path, tmp_path: Path, vcr_cassette
):
    vcr_cassette.allow_playback_repeats = False
    output_vcf = tmp_path / "test_vcf_output_grch38_noattrs.vcf.gz"
    output_vrs_pkl = tmp_path / "test_vcf_pkl_grch38_noattrs.pkl"
    expected_vcf_no_vrs_attrs = (
        TEST_DATA_DIR / "test_vcf_expected_output_no_vrs_attrs.vcf"
    )

    # Test GRCh38 assembly, which was used for input_vcf and no vrs attributes
    vcf_annotator.annotate(input_vcf, output_vcf, output_pkl_path=output_vrs_pkl)
    compare_vcfs(output_vcf, expected_vcf_no_vrs_attrs)
    assert output_vrs_pkl.exists()
    if vcr_cassette.write_protected:
        assert vcr_cassette.all_played


@pytest.mark.vcr
def test_annotate_vcf_grch38_attrs(
    vcf_annotator: VcfAnnotator, input_vcf: Path, tmp_path: Path, vcr_cassette
):
    vcr_cassette.allow_playback_repeats = False
    output_vcf = tmp_path / "test_vcf_output_grch38_attrs.vcf.gz"
    output_vrs_pkl = tmp_path / "test_vcf_pkl_grch38_attrs.pkl"
    expected_vcf = TEST_DATA_DIR / "test_vcf_expected_output.vcf"

    # Test GRCh38 assembly, which was used for input_vcf and vrs attributes
    vcf_annotator.annotate(
        input_vcf, output_vcf, vrs_attributes=True, output_pkl_path=output_vrs_pkl
    )
    compare_vcfs(output_vcf, expected_vcf)
    assert output_vrs_pkl.exists()
    if vcr_cassette.write_protected:
        assert vcr_cassette.all_played


@pytest.mark.vcr
def test_annotate_vcf_grch38_attrs_altsonly(
    vcf_annotator: VcfAnnotator, input_vcf: Path, tmp_path: Path, vcr_cassette
):
    vcr_cassette.allow_playback_repeats = False
    output_vcf = tmp_path / "test_vcf_output_grch38_attrs_altsonly.vcf.gz"
    output_vrs_pkl = tmp_path / "test_vcf_pkl_grch38_attrs_altsonly.pkl"
    expected_altsonly_vcf = TEST_DATA_DIR / "test_vcf_expected_altsonly_output.vcf"

    # Test GRCh38 assembly with VRS computed for ALTs only, which was used for input_vcf and vrs attributes
    vcf_annotator.annotate(
        input_vcf,
        output_vcf,
        vrs_attributes=True,
        compute_for_ref=False,
        output_pkl_path=output_vrs_pkl,
    )
    compare_vcfs(output_vcf, expected_altsonly_vcf)
    assert output_vrs_pkl.exists()
    if vcr_cassette.write_protected:
        assert vcr_cassette.all_played


@pytest.mark.vcr
def test_annotate_vcf_grch37_attrs(
    vcf_annotator: VcfAnnotator, input_vcf: Path, tmp_path: Path, vcr_cassette
):
    vcr_cassette.allow_playback_repeats = False
    output_vcf = tmp_path / "test_vcf_output_grch37_attrs.vcf.gz"
    output_vrs_pkl = tmp_path / "test_vcf_pkl_grch37_attrs.pkl"
    expected_vcf = TEST_DATA_DIR / "test_vcf_expected_output.vcf"

    # Test GRCh37 assembly, which was not used for input_vcf
    vcf_annotator.annotate(
        input_vcf,
        output_vcf,
        vrs_attributes=True,
        assembly="GRCh37",
        output_pkl_path=output_vrs_pkl,
    )
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    with expected_vcf.open() as expected_output:
        expected_output_lines = expected_output.readlines()
    assert out_vcf_lines != expected_output_lines
    assert output_vrs_pkl.exists()
    if vcr_cassette.write_protected:
        assert vcr_cassette.all_played


@pytest.mark.vcr
def test_annotate_vcf_pickle_only(
    vcf_annotator: VcfAnnotator, input_vcf: Path, tmp_path: Path, vcr_cassette
):
    vcr_cassette.allow_playback_repeats = False
    output_vcf = tmp_path / "test_vcf_output_pickle_only.vcf.gz"
    output_vrs_pkl = tmp_path / "test_vcf_pkl_pickle_only.pkl"

    # Test only pickle output
    vcf_annotator.annotate(
        input_vcf, output_pkl_path=output_vrs_pkl, vrs_attributes=True
    )
    assert output_vrs_pkl.exists()
    assert not output_vcf.exists()
    if vcr_cassette.write_protected:
        assert vcr_cassette.all_played


@pytest.mark.vcr
def test_annotate_vcf_vcf_only(
    vcf_annotator: VcfAnnotator, input_vcf: Path, tmp_path: Path, vcr_cassette
):
    vcr_cassette.allow_playback_repeats = False
    output_vcf = tmp_path / "test_vcf_output_vcf_only.vcf.gz"
    output_vrs_pkl = tmp_path / "test_vcf_pkl_vcf_only.pkl"
    expected_vcf = TEST_DATA_DIR / "test_vcf_expected_output.vcf"

    # Test only VCF output
    vcf_annotator.annotate(input_vcf, output_vcf_path=output_vcf, vrs_attributes=True)
    compare_vcfs(output_vcf, expected_vcf)
    if vcr_cassette.write_protected:
        assert vcr_cassette.all_played
    assert not Path(output_vrs_pkl).exists()


def test_annotate_vcf_input_validation(vcf_annotator: VcfAnnotator, input_vcf: Path):
    with pytest.raises(
        VcfAnnotatorError,
        match="No VCF, PKL, or NDJSON output path provided -- must pass at least one of `output_vcf_path`, `output_pkl_path`, `output_ndjson_path` to annotate().",
    ):
        vcf_annotator.annotate(input_vcf)


@pytest.mark.vcr
def test_get_vrs_object_invalid_input(vcf_annotator: VcfAnnotator, caplog):
    """Test that _get_vrs_object method works as expected with invalid input"""
    # some tests below are checking for debug logging statements
    caplog.set_level(logging.DEBUG)

    # No CHROM
    vcf_annotator._get_vrs_object(".-140753336-A-T", [], {}, "GRCh38")
    assert "KeyError when getting refget accession: GRCh38:." in caplog.text

    # No POS
    vcf_annotator._get_vrs_object("7-.-A-T", [], {}, "GRCh38")
    assert "None was returned when translating 7-.-A-T from gnomad" in caplog.text

    # No REF
    vcf_annotator._get_vrs_object("7-140753336-.-T", [], {}, "GRCh38")
    assert (
        "None was returned when translating 7-140753336-.-T from gnomad" in caplog.text
    )

    # No ALT
    vcf_annotator._get_vrs_object("7-140753336-A-.", [], {}, "GRCh38")
    assert (
        "None was returned when translating 7-140753336-A-. from gnomad" in caplog.text
    )

    # Invalid ref, but not requiring validation checks so no error is raised
    vcf_annotator._get_vrs_object(
        "7-140753336-G-T", [], {}, "GRCh38", require_validation=False
    )
    assert "" in caplog.text

    # Invalid ref, but requiring validation checks so an error is raised
    invalid_ref_seq_msg = "Reference mismatch at GRCh38:7 position 140753335-140753336 (input gave 'C' but correct ref is 'A')"
    with pytest.raises(DataProxyValidationError, match=re.escape(invalid_ref_seq_msg)):
        vcf_annotator._get_vrs_object(
            "7-140753336-C-T", [], {}, "GRCh38", require_validation=True
        )
    assert invalid_ref_seq_msg in caplog.text


@pytest.mark.vcr
def test_annotate_vcf_rle(vcf_annotator: VcfAnnotator, vcr_cassette):
    """Test VCF annotation with ReferenceLengthExpression states.

    Tests two variants that should produce ReferenceLengthExpression states:
    1. Deletion (1:100210778 AA>A): Expected RLE with length=1, repeatSubunitLength=1
    2. Duplication (1:102995989 CTTT>CTTTCTTT): Expected RLE with length=8, repeatSubunitLength=4
    """
    vcr_cassette.allow_playback_repeats = False
    input_vcf = TEST_DATA_DIR / "test_rle.vcf"
    output_vcf = TEST_DATA_DIR / "test_rle_output.vcf"
    output_vrs_pkl = TEST_DATA_DIR / "test_rle_output.pkl"

    # Annotate the VCF with VRS attributes enabled
    vcf_annotator.annotate(
        input_vcf, output_vcf, vrs_attributes=True, output_pkl_path=output_vrs_pkl
    )

    # Read the output VCF and verify RLE fields are present
    with pysam.VariantFile(str(output_vcf)) as vcf:
        # Verify the RLE-specific header fields were added
        assert "VRS_Lengths" in vcf.header.info
        assert "VRS_RepeatSubunitLengths" in vcf.header.info

        variants = list(vcf)
        assert len(variants) == 2

        # Test variant 1: Deletion (AA>A)
        # Expected: length=1, repeatSubunitLength=1
        deletion_variant = variants[0]
        assert deletion_variant.chrom == "1"
        assert deletion_variant.pos == 100210778
        assert deletion_variant.ref == "AA"
        assert deletion_variant.alts == ("A",)

        # Check VRS attributes for deletion
        assert "VRS_Allele_IDs" in deletion_variant.info
        assert "VRS_Starts" in deletion_variant.info
        assert "VRS_Ends" in deletion_variant.info
        assert "VRS_States" in deletion_variant.info
        assert "VRS_Lengths" in deletion_variant.info
        assert "VRS_RepeatSubunitLengths" in deletion_variant.info

        # Expected values for deletion RLE
        # REF: AA uses RLE with length=2, repeatSubunitLength=2
        # ALT: A uses RLE with length=1, repeatSubunitLength=1
        vrs_lengths = deletion_variant.info["VRS_Lengths"]
        vrs_repeat_lengths = deletion_variant.info["VRS_RepeatSubunitLengths"]
        assert len(vrs_lengths) == 2  # REF and ALT
        assert len(vrs_repeat_lengths) == 2
        assert vrs_lengths == (2, 1)  # REF: AA (length 2), ALT: A (length 1)
        assert vrs_repeat_lengths == (
            2,
            1,
        )  # REF: AA as repeat unit of length 2, ALT: A as repeat unit of length 1

        # Test variant 2: Duplication (CTTT>CTTTCTTT)
        # Expected: length=8, repeatSubunitLength=4
        duplication_variant = variants[1]
        assert duplication_variant.chrom == "1"
        assert duplication_variant.pos == 102995989
        assert duplication_variant.ref == "CTTT"
        assert duplication_variant.alts == ("CTTTCTTT",)

        # Check VRS attributes for duplication
        assert "VRS_Allele_IDs" in duplication_variant.info
        assert "VRS_Lengths" in duplication_variant.info
        assert "VRS_RepeatSubunitLengths" in duplication_variant.info

        # Expected values for duplication RLE
        # REF: CTTT uses RLE with length=4, repeatSubunitLength=4
        # ALT: CTTTCTTT uses RLE with length=8, repeatSubunitLength=4
        vrs_lengths = duplication_variant.info["VRS_Lengths"]
        vrs_repeat_lengths = duplication_variant.info["VRS_RepeatSubunitLengths"]
        assert len(vrs_lengths) == 2  # REF and ALT
        assert len(vrs_repeat_lengths) == 2
        assert vrs_lengths == (4, 8)  # REF: CTTT (length 4), ALT: CTTTCTTT (length 8)
        assert vrs_repeat_lengths == (4, 4)  # Both are 4-base repeats

    assert output_vrs_pkl.exists()
    assert vcr_cassette.all_played
