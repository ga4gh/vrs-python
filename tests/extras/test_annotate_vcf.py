"""Ensure proper functionality of VCFAnnotator"""

import gzip
import logging
import os
import re
from pathlib import Path

import pytest

from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy

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


@pytest.mark.vcr
def test_annotate_vcf_grch38_noattrs(
    vcf_annotator: VcfAnnotator, input_vcf: Path, tmp_path: Path, vcr_cassette
):
    vcr_cassette.allow_playback_repeats = False
    output_vcf = tmp_path / "test_vcf_output_grch38_noattrs.vcf.gz"
    output_vrs_pkl = tmp_path / "test_vcf_pkl_grch38_noattrs.pkl"
    expected_vcf_no_vrs_attrs = (
        TEST_DATA_DIR / "test_vcf_expected_output_no_vrs_attrs.vcf.gz"
    )

    # Test GRCh38 assembly, which was used for input_vcf and no vrs attributes
    vcf_annotator.annotate(input_vcf, output_vcf, output_vrs_pkl)
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    with gzip.open(expected_vcf_no_vrs_attrs, "rt") as expected_output:
        expected_output_lines = expected_output.readlines()
    for actual_line, expected_line in zip(
        out_vcf_lines, expected_output_lines, strict=False
    ):
        assert actual_line == expected_line
    assert output_vrs_pkl.exists()
    assert vcr_cassette.all_played


@pytest.mark.vcr
def test_annotate_vcf_grch38_attrs(
    vcf_annotator: VcfAnnotator, input_vcf: Path, tmp_path: Path, vcr_cassette
):
    vcr_cassette.allow_playback_repeats = False
    output_vcf = tmp_path / "test_vcf_output_grch38_attrs.vcf.gz"
    output_vrs_pkl = tmp_path / "test_vcf_pkl_grch38_attrs.pkl"
    expected_vcf = TEST_DATA_DIR / "test_vcf_expected_output.vcf.gz"

    # Test GRCh38 assembly, which was used for input_vcf and vrs attributes
    vcf_annotator.annotate(input_vcf, output_vcf, output_vrs_pkl, vrs_attributes=True)
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    with gzip.open(expected_vcf, "rt") as expected_output:
        expected_output_lines = expected_output.readlines()
    for actual_line, expected_line in zip(
        out_vcf_lines, expected_output_lines, strict=False
    ):
        assert actual_line == expected_line
    assert output_vrs_pkl.exists()
    assert vcr_cassette.all_played


@pytest.mark.vcr
def test_annotate_vcf_grch38_attrs_altsonly(
    vcf_annotator: VcfAnnotator, input_vcf: Path, tmp_path: Path, vcr_cassette
):
    vcr_cassette.allow_playback_repeats = False
    output_vcf = tmp_path / "test_vcf_output_grch38_attrs_altsonly.vcf.gz"
    output_vrs_pkl = tmp_path / "test_vcf_pkl_grch38_attrs_altsonly.pkl"
    expected_altsonly_vcf = TEST_DATA_DIR / "test_vcf_expected_altsonly_output.vcf.gz"

    # Test GRCh38 assembly with VRS computed for ALTs only, which was used for input_vcf and vrs attributes
    vcf_annotator.annotate(
        input_vcf,
        output_vcf,
        output_vrs_pkl,
        vrs_attributes=True,
        compute_for_ref=False,
    )
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    with gzip.open(expected_altsonly_vcf, "rt") as expected_output:
        expected_output_lines = expected_output.readlines()
    for actual_line, expected_line in zip(
        out_vcf_lines, expected_output_lines, strict=False
    ):
        assert actual_line == expected_line
    assert output_vrs_pkl.exists()
    assert vcr_cassette.all_played


@pytest.mark.vcr
def test_annotate_vcf_grch37_attrs(
    vcf_annotator: VcfAnnotator, input_vcf: Path, tmp_path: Path, vcr_cassette
):
    vcr_cassette.allow_playback_repeats = False
    output_vcf = tmp_path / "test_vcf_output_grch37_attrs.vcf.gz"
    output_vrs_pkl = tmp_path / "test_vcf_pkl_grch37_attrs.pkl"
    expected_vcf = TEST_DATA_DIR / "test_vcf_expected_output.vcf.gz"

    # Test GRCh37 assembly, which was not used for input_vcf
    vcf_annotator.annotate(
        input_vcf, output_vcf, output_vrs_pkl, vrs_attributes=True, assembly="GRCh37"
    )
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    with gzip.open(expected_vcf, "rt") as expected_output:
        expected_output_lines = expected_output.readlines()
    assert out_vcf_lines != expected_output_lines
    assert output_vrs_pkl.exists()
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
    assert vcr_cassette.all_played


@pytest.mark.vcr
def test_annotate_vcf_vcf_only(
    vcf_annotator: VcfAnnotator, input_vcf: Path, tmp_path: Path, vcr_cassette
):
    vcr_cassette.allow_playback_repeats = False
    output_vcf = tmp_path / "test_vcf_output_vcf_only.vcf.gz"
    output_vrs_pkl = tmp_path / "test_vcf_pkl_vcf_only.pkl"
    expected_vcf = TEST_DATA_DIR / "test_vcf_expected_output.vcf.gz"

    # Test only VCF output
    vcf_annotator.annotate(input_vcf, output_vcf_path=output_vcf, vrs_attributes=True)
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    with gzip.open(expected_vcf, "rt") as expected_output:
        expected_output_lines = expected_output.readlines()
    assert out_vcf_lines == expected_output_lines
    assert vcr_cassette.all_played
    assert not Path(output_vrs_pkl).exists()


def test_annotate_vcf_input_validation(vcf_annotator: VcfAnnotator, input_vcf: Path):
    with pytest.raises(VCFAnnotatorError) as e:
        vcf_annotator.annotate(input_vcf)
    assert str(e.value) == "Must provide one of: `vcf_out` or `vrs_pickle_out`"


@pytest.mark.vcr
def test_get_vrs_object_invalid_input(vcf_annotator: VcfAnnotator, caplog):
    """Test that _get_vrs_object method works as expected with invalid input"""
    # some tests below are checking for debug logging statements
    caplog.set_level(logging.DEBUG)

    # No CHROM
    vcf_annotator._get_vrs_object(".-140753336-A-T", {}, [], "GRCh38")
    assert "KeyError when getting refget accession: GRCh38:." in caplog.text

    # No POS
    vcf_annotator._get_vrs_object("7-.-A-T", {}, [], "GRCh38")
    assert "None was returned when translating 7-.-A-T from gnomad" in caplog.text

    # No REF
    vcf_annotator._get_vrs_object("7-140753336-.-T", {}, [], "GRCh38")
    assert (
        "None was returned when translating 7-140753336-.-T from gnomad" in caplog.text
    )

    # No ALT
    vcf_annotator._get_vrs_object("7-140753336-A-.", {}, [], "GRCh38")
    assert (
        "None was returned when translating 7-140753336-A-. from gnomad" in caplog.text
    )

    # Invalid ref, but not requiring validation checks so no error is raised
    vcf_annotator._get_vrs_object(
        "7-140753336-G-T", {}, [], "GRCh38", require_validation=False
    )
    assert "" in caplog.text

    # Invalid ref, but requiring validation checks so an error is raised
    invalid_ref_seq_msg = "Expected reference sequence C on GRCh38:7 at positions (140753335, 140753336) but found A"
    with pytest.raises(DataProxyValidationError, match=re.escape(invalid_ref_seq_msg)):
        vcf_annotator._get_vrs_object(
            "7-140753336-C-T", {}, [], "GRCh38", require_validation=True
        )
    assert invalid_ref_seq_msg in caplog.text
