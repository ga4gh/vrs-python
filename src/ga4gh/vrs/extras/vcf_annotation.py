"""
Annotate VCF files with VRS

Input Format: VCF
Output Format: VCF

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
    Annotates an input VCF file with VRS allele ids & creates a 
    pickle file containing the vrs object information

    """

    def __init__(self, tlr) -> None:
        self.tlr = tlr

    def annotate(self, inputfile, outputfile, vrsfile):
        INFO_FIELD_ID = "VRS_Allele"
        vrs_data = {}
        vcf_in = pysam.VariantFile(filename=inputfile)
        vcf_in.header.info.add(INFO_FIELD_ID, "1", "String", "vrs")
        vcf_out = pysam.VariantFile(outputfile, "w", header=vcf_in.header)
        vrs_out = open(vrsfile, "wb")    #For sending VRS data to the pickle file

        for record in vcf_in:
            ld = self._record_digests(record, vrs_data)
            record.info[INFO_FIELD_ID] = ",".join(ld)
            vcf_out.write(record)

        pickle.dump(vrs_data, vrs_out)

        vrs_out.close()
        vcf_in.close()

    def _record_digests(self, record, vrs_data):
        gnomad_loc = f"{record.chrom}-{record.pos}"
        alts = record.alts if record.alts else []
        data = f"{record.chrom}\t{record.pos}\t{record.id}\t{record.ref}\t{record.alts}"
        # payloads like ['20:14369:1', '20:14369:1:G', '20:14369:1:A']
        reference_allele = f"{gnomad_loc}-{record.ref}-{record.ref}"
        vrs_ref_object = self.tlr.translate_from(reference_allele, "gnomad")
        vrs_ref_object_id = vrs_ref_object._id._value
        alleles = [f"{gnomad_loc}-{record.ref}-{a}" for a in [*alts]]    #using gnomad format
        vrs_allele = [vrs_ref_object_id]
        for allele in alleles:
            if "*" in allele:
                vrs_allele.append("")
            else:
                vrs_object = self.tlr.translate_from(allele, "gnomad")
                vrs_allele.append(vrs_object._id._value)
                vrs_data[data] = str(vrs_object)

        return vrs_allele


def parse_args(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("VCF_IN")
    ap.add_argument("--out", "-o", default="-")
    ap.add_argument("--vrs-file", default="-")
    opts = ap.parse_args(argv)
    return opts


if __name__ == "__main__":
    start_time = time.time()

    opts = parse_args(sys.argv[1:])
    print(f"These are the options that you have selected: {opts}\n")
    data_proxy = SeqRepoDataProxy(SeqRepo("/usr/local/share/seqrepo/latest"))
    tlr = Translator(data_proxy)
    vcf_annotator = VCFAnnotator(tlr)
    vcf_annotator.annotate(opts.VCF_IN, opts.out, opts.vrs_file)

    end_time = time.time()
    total_time = (float(end_time) - float(start_time))
    total_time_minutes = (total_time / 60)
    print(f"This program took {total_time} seconds to run.")
    print(f"This program took {total_time_minutes} minutes to run.")
