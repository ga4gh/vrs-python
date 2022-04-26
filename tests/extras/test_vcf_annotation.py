"""
Purpose: Ensure proper functionality of VCF_Annotation.py

Input Format: VCF
Output Format: VCF

Function 1: Test proper vrs_allele_id creation
    - 
Function 2: Test addition of data to VCF file

"""
import pytest
from biocommons.seqrepo import SeqRepo
import pysam
from ga4gh.vrs.extras.translator import Translator
from ga4gh.vrs.extras.vcf_annotation import VCFAnnotator


@pytest.mark.vcr
def test_annotate_vcf(tlr):
    annotate_vcf("tests/extras/test_vcf_annotation.py")


def test_record_digests():
    resp = record_digests
