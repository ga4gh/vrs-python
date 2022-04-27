"""
Purpose: Ensure proper functionality of VCF_Annotation.py

Input Format: VCF
Output Format: VCF

Function 1: Test proper vrs_allele_id creation
    - 
Function 2: Test addition of data to VCF file

"""
import gzip
import os

import pytest

from ga4gh.vrs.extras.vcf_annotation import VCFAnnotator
from ga4gh.vrs.extras.translator import Translator


@pytest.fixture(scope="module")
def tlr_local(dataproxy):
    return Translator(data_proxy=dataproxy)


@pytest.fixture(scope="module")
def vcf_annotator(tlr_local):
    return VCFAnnotator(tlr_local)


@pytest.mark.vcr
def test_annotate_vcf(vcf_annotator):
    input_vcf = "tests/extras/data/test_vcf.vcf"
    output_vcf = "tests/extras/data/test_vcf_out.vcf.gz"
    output_vrs_pkl = "tests/extras/data/test_vcf_pkl.pkl"
    expected_vcf = "tests/extras/data/test_expected_output_vcf.vcf"
    vcf_annotator.annotate(input_vcf, output_vcf, output_vrs_pkl)
    with gzip.open(output_vcf, "rb") as out_vcf:
        #print(out_vcf.readlines())
        out_vcf_lines = out_vcf.readlines()
    with open(expected_vcf, "r") as expected_output:
        expected_output_lines = expected_output.readlines()
        assert out_vcf_lines == expected_output_lines

    assert os.path.exists(output_vrs_pkl)

    #os.remove(out_vcf)
    os.remove(output_vrs_pkl)
