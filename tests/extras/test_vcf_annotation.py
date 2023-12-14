"""Ensure proper functionality of VCFAnnotator"""
import gzip
import os

import pytest

from ga4gh.vrs.extras.vcf_annotation import VCFAnnotator, VCFAnnotatorException


@pytest.fixture(scope="module")
def vcf_annotator():
    return VCFAnnotator("rest")

@pytest.mark.vcr
def test_annotate_vcf(vcf_annotator):
    TEST_DATA_DIR = "tests/extras/data"

    input_vcf = f"{TEST_DATA_DIR}/test_vcf_input.vcf"
    output_vcf = f"{TEST_DATA_DIR}/test_vcf_output.vcf.gz"
    output_vrs_pkl = f"{TEST_DATA_DIR}/test_vcf_pkl.pkl"
    expected_vcf = f"{TEST_DATA_DIR}/test_vcf_expected_output.vcf.gz"
    expected_altsonly_vcf = f"{TEST_DATA_DIR}/test_vcf_expected_altsonly_output.vcf.gz"
    expected_vcf_no_vrs_attrs = f"{TEST_DATA_DIR}/test_vcf_expected_output_no_vrs_attrs.vcf.gz"

    # Test GRCh38 assembly, which was used for input_vcf and no vrs attributes
    vcf_annotator.annotate(input_vcf, output_vcf, output_vrs_pkl)
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    with gzip.open(expected_vcf_no_vrs_attrs, "rt") as expected_output:
        expected_output_lines = expected_output.readlines()
    assert out_vcf_lines == expected_output_lines
    assert os.path.exists(output_vrs_pkl)
    os.remove(output_vcf)
    os.remove(output_vrs_pkl)

    # Test GRCh38 assembly, which was used for input_vcf and vrs attributes
    vcf_annotator.annotate(input_vcf, output_vcf, output_vrs_pkl, vrs_attributes=True)
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    with gzip.open(expected_vcf, "rt") as expected_output:
        expected_output_lines = expected_output.readlines()
    assert out_vcf_lines == expected_output_lines
    assert os.path.exists(output_vrs_pkl)
    os.remove(output_vcf)
    os.remove(output_vrs_pkl)

    # Test GRCh38 assembly with VRS computed for ALTs only, which was used for input_vcf and vrs attributes
    vcf_annotator.annotate(input_vcf, output_vcf, output_vrs_pkl, vrs_attributes=True, compute_for_ref=False)
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    with gzip.open(expected_altsonly_vcf, "rt") as expected_output:
        expected_output_lines = expected_output.readlines()
    assert out_vcf_lines == expected_output_lines
    assert os.path.exists(output_vrs_pkl)
    os.remove(output_vcf)
    os.remove(output_vrs_pkl)

    # Test GRCh37 assembly, which was not used for input_vcf
    vcf_annotator.annotate(input_vcf, output_vcf, output_vrs_pkl, vrs_attributes=True, assembly="GRCh37")
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    with gzip.open(expected_vcf, "rt") as expected_output:
        expected_output_lines = expected_output.readlines()
    assert out_vcf_lines != expected_output_lines
    assert os.path.exists(output_vrs_pkl)
    os.remove(output_vcf)
    os.remove(output_vrs_pkl)

    # Test only pickle output
    vcf_annotator.annotate(input_vcf, vrs_pickle_out=output_vrs_pkl, vrs_attributes=True)
    assert os.path.exists(output_vrs_pkl)
    assert not os.path.exists(output_vcf)
    os.remove(output_vrs_pkl)

    # Test only VCF output
    vcf_annotator.annotate(input_vcf, vcf_out=output_vcf, vrs_attributes=True)
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    assert out_vcf_lines == expected_output_lines
    os.remove(output_vcf)
    assert not os.path.exists(output_vrs_pkl)

    with pytest.raises(VCFAnnotatorException) as e:
        vcf_annotator.annotate(input_vcf)
    assert str(e.value) == "Must provide one of: `vcf_out` or `vrs_pickle_out`"

@pytest.mark.vcr
def test_get_vrs_object_invalid_input(vcf_annotator, caplog):
    """Test that _get_vrs_object method works as expected with invalid input"""
    # No CHROM
    vcf_annotator._get_vrs_object(".-140753336-A-T", {}, [], "GRCh38")
    assert "KeyError when getting refget accession: GRCh38:." in caplog.text

    # No POS
    vcf_annotator._get_vrs_object("7-.-A-T", {}, [], "GRCh38")
    assert "None was returned when translating 7-.-A-T from gnomad" in caplog.text

    # No REF
    vcf_annotator._get_vrs_object("7-140753336-.-T", {}, [], "GRCh38")
    assert "None was returned when translating 7-140753336-.-T from gnomad" in caplog.text

    # No ALT
    vcf_annotator._get_vrs_object("7-140753336-A-.", {}, [], "GRCh38")
    assert "None was returned when translating 7-140753336-A-. from gnomad" in caplog.text
