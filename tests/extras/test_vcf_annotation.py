"""Ensure proper functionality of VCFAnnotator"""
import gzip
import os
import tempfile

import pytest

from ga4gh.vrs.extras.vcf_annotation import VCFAnnotator
from ga4gh.vrs.extras.translator import Translator


@pytest.fixture(scope="module")
def vcf_annotator():
    root_dir = os.environ.get("SEQREPO_ROOT_DIR", "/usr/local/share/seqrepo/latest")
    return VCFAnnotator("local", seqrepo_root_dir=root_dir)

@pytest.mark.vcr
def test_annotate_vcf(vcf_annotator):
    TEST_DATA_DIR = "tests/extras/data"

    input_vcf = f"{TEST_DATA_DIR}/test_vcf_input.vcf"
    output_vcf = f"{tempfile.gettempdir()}/test_vcf_output.vcf.gz"
    output_vrs_pkl = f"{tempfile.gettempdir()}/test_vcf_pkl.pkl"
    expected_vcf = f"{TEST_DATA_DIR}/test_vcf_expected_output.vcf.gz"

    # Test GRCh38 assembly, which was used for input_vcf
    vcf_annotator.annotate(input_vcf, output_vcf, output_vrs_pkl)
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    with gzip.open(expected_vcf, "rt") as expected_output:
        expected_output_lines = expected_output.readlines()
    assert out_vcf_lines == expected_output_lines
    assert os.path.exists(output_vrs_pkl)
    os.remove(output_vcf)
    os.remove(output_vrs_pkl)

    # Test GRCh37 assembly, which was not used for input_vcf
    vcf_annotator.annotate(input_vcf, output_vcf, output_vrs_pkl, "GRCh37")
    with gzip.open(output_vcf, "rt") as out_vcf:
        out_vcf_lines = out_vcf.readlines()
    assert out_vcf_lines != expected_output_lines
    assert os.path.exists(output_vrs_pkl)
    os.remove(output_vcf)
    os.remove(output_vrs_pkl)
