"""
Annotate VCF files with VRS

Input Format: VCF
Output Format: VCF

The user should pass arguments for the VCF input, VCF output, &
the vrs object file name.

ex. python3 src/ga4gh/vrs/extras/vcf_annotation.py input.vcf.gz --out
./output.vcf.gz --vrs-file ./vrs_objects.pkl
"""

import argparse
import sys
import pickle
import time

from biocommons.seqrepo import SeqRepo
import pysam

from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator


class VCFAnnotator:
    """
    This class provides utility for annotating VCF's with VRS allele id's.

    VCF's are read using pysam and stored as pysam objects.
    Alleles are translated into vrs allele id's using VRS-Python Translator.

    """

    def __init__(self, tlr) -> None:
        """
        param: Translator tlr Valid translator object with a specified data proxy
        """
        self.tlr = tlr

    def annotate(self, inputfile, outputfile, vrsfile):
        """
        Annotates an input VCF file with VRS allele ids & creates a
        pickle file containing the vrs object information.
        param: str inputfile The path and filename for the input VCF file
        param: str outputfile The path and filename for the output VCF file
        param: str vrsfile The path and filename for the output VRS object file
        """
        INFO_FIELD_ID = "VRS_Allele"
        vrs_data = {}
        vcf_in = pysam.VariantFile(filename=inputfile)
        vcf_in.header.info.add(INFO_FIELD_ID, "1", "String", "vrs")
        vcf_out = pysam.VariantFile(outputfile, "w", header=vcf_in.header)
        vrs_out = open(vrsfile, "wb")    # For sending VRS data to the pickle file

        for record in vcf_in:
            ld = self._record_digests(record, vrs_data)
            record.info[INFO_FIELD_ID] = ",".join(ld)
            vcf_out.write(record)

        pickle.dump(vrs_data, vrs_out)

        vrs_out.close()
        vcf_in.close()
        vcf_out.close()

    def _record_digests(self, record, vrs_data):
        """
        Mutate vrs_data with vrs object information and returning a list of vrs allele ids
        param: pysam.VariantRecord record A row in the vcf file
        param: dict vrs_data Dictionary containing the VRS object information for the VCF
        return: list vrs_allele_ids List containing the vrs allele id information
        """
        gnomad_loc = f"{record.chrom}-{record.pos}"
        alts = record.alts if record.alts else []
        data = f"{record.chrom}\t{record.pos}\t{record.id}\t{record.ref}\t{record.alts}"
        # payloads like ['20:14369:1', '20:14369:1:G', '20:14369:1:A']
        reference_allele = f"{gnomad_loc}-{record.ref}-{record.ref}"
        vrs_ref_object = self.tlr.translate_from(reference_allele, "gnomad")
        vrs_data[reference_allele] = str(vrs_ref_object.as_dict())
        alleles = [f"{gnomad_loc}-{record.ref}-{a}" for a in [*alts]]    # using gnomad format
        vrs_allele_ids = [vrs_ref_object.id._value]
        for allele in alleles:
            if "*" in allele:
                vrs_allele_ids.append("")
            else:
                vrs_object = self.tlr.translate_from(allele, "gnomad")
                vrs_allele_ids.append(vrs_object.id._value)
                vrs_data[data] = str(vrs_object.as_dict())

        return vrs_allele_ids


def parse_args(argv):
    """
    Parses arguments passed in by the user
    param: list[str] argv Arguments passed by the user to specify file locations and names
    return: argparse.Namespace Returns the options passed by the user to be assigned to proper variables
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("VCF_IN")
    ap.add_argument("--out", "-o", default="-")
    ap.add_argument("--vrs-file", default="-")
    opts = ap.parse_args(argv)
    return opts


if __name__ == "__main__":
    start_time = time.time()

    options = parse_args(sys.argv[1:])
    print(f"These are the options that you have selected: {options}\n")
    data_proxy = SeqRepoDataProxy(SeqRepo("/usr/local/share/seqrepo/latest"))
    tlr = Translator(data_proxy)
    vcf_annotator = VCFAnnotator(tlr)
    vcf_annotator.annotate(options.VCF_IN, options.out, options.vrs_file)

    end_time = time.time()
    total_time = (float(end_time) - float(start_time))
    total_time_minutes = (total_time / 60)
    print(f"This program took {total_time} seconds to run.")
    print(f"This program took {total_time_minutes} minutes to run.")
