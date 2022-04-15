"""
Annotate VCF files with VRS

Input Format: VCF
Output Format: VCF

"""

import argparse
import hashlib
import logging
import sys
import pickle
import time
import coloredlogs
from biocommons.seqrepo import SeqRepo
from pysam import VariantFile
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator

dp = SeqRepoDataProxy(SeqRepo("/usr/local/share/seqrepo/latest"))
tlr = Translator(data_proxy=dp)

_logger = logging.getLogger()

INFO_FIELD_ID = "VRS_Allele"
vrs_data = {}


def parse_args(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("VCF_IN")
    ap.add_argument("--out", "-o", default="-")
    ap.add_argument("--failure-probability", "-f", type=float, default=0)
    ap.add_argument("--vrs_file", default="-")
    opts = ap.parse_args(argv)
    return opts


def record_digests(record):

    def digest(d):
        return hashlib.sha512(d.encode("ascii")).hexdigest()[:8]

    #spdi_loc = f"{record.chrom}:{record.pos-1}:{len(record.ref)}"
    gnomad_loc = f"{record.chrom}-{record.pos}"
    alts = record.alts if record.alts else []
    data = f"{record.chrom}\t{record.pos}\t{record.id}\t{record.ref}\t{record.alts}"
    # payloads like ['20:14369:1', '20:14369:1:G', '20:14369:1:A']
    reference_allele = f"{gnomad_loc}-{record.ref}-{record.ref}"
    vrs_ref_object = tlr.translate_from(reference_allele, "gnomad")
    vrs_ref_object_id = vrs_ref_object._id._value
    alleles = [f"{gnomad_loc}-{record.ref}-{a}" for a in [*alts]]    #using gnomad format
    vrs_allele = []
    vrs_allele.append(vrs_ref_object_id)
    for allele in alleles:
        if "*" in allele:
            vrs_allele.append("intentional_blank")
        else:
            vrs_object = tlr.translate_from(allele, "gnomad")
            vrs_allele.append(vrs_object._id._value)
            vrs_data[data] = str(vrs_object)

    return vrs_allele


if __name__ == "__main__":
    start_time = time.time()
    coloredlogs.install(level="INFO")

    opts = parse_args(sys.argv[1:])
    print(f"These are the options that you have selected: {opts}")

    vcf_in = VariantFile(filename=opts.VCF_IN)
    vcf_in.header.info.add(INFO_FIELD_ID, "1", "String", "vrs")
    vcf_out = VariantFile(opts.out, "w", header=vcf_in.header)
    vrs_out = open(opts.vrs_file, "wb")    #For sending VRS data to the pickle file

    for record in vcf_in:
        ld = record_digests(record)
        record.info[INFO_FIELD_ID] = ",".join(ld)
        vcf_out.write(record)

    #print(vrs_data)
    pickle.dump(vrs_data, vrs_out)

    vrs_out.close()
    vcf_in.close()
    end_time = time.time()
    total_time = (float(end_time) - float(start_time))
    total_time_minutes = (total_time / 60)
    print(f"This program took {total_time} seconds to run.")
    print(f"This program took {total_time_minutes} minutes to run.")
