"""Module containing tools for annotating VCFs with VRS

Example of how to run from root of vrs-python directory:
python3 -m src.ga4gh.vrs.extras.vcf_annotation --vcf_in input.vcf.gz \
    --vcf_out output.vcf.gz --vrs_file vrs_objects.pkl
"""
import pickle
from enum import Enum
from typing import Dict, List, Optional
from timeit import default_timer as timer

import click
from biocommons.seqrepo import SeqRepo
import pysam

from ga4gh.vrs.dataproxy import SeqRepoDataProxy, SeqRepoRESTDataProxy
from ga4gh.vrs.extras.translator import Translator


class SeqRepoProxyType(str, Enum):
    """Define constraints for SeqRepo Data Proxy types"""

    LOCAL = "local"
    REST = "rest"


@click.command()
@click.option(
    "--vcf_in",
    required=True,
    type=str,
    help="The path for the input VCF file to annotate"
)
@click.option(
    "--vcf_out",
    required=True,
    type=str,
    help="The path for the output VCF file"
)
@click.option(
    "--vrs_file",
    required=True,
    type=str,
    help="The path for the output VCF pickle file"
)
@click.option(
    "--seqrepo_dp_type",
    required=False,
    default=SeqRepoProxyType.LOCAL,
    type=click.Choice([v.value for v in SeqRepoProxyType.__members__.values()],
                      case_sensitive=True),
    help="The type of the SeqRepo Data Proxy to use",
    show_default=True,
    show_choices=True
)
@click.option(
    "--seqrepo_root_dir",
    required=False,
    default="/usr/local/share/seqrepo/latest",
    help="The root directory for local SeqRepo instance",
    show_default=True
)
@click.option(
    "--seqrepo_base_url",
    required=False,
    default="http://localhost:5000/seqrepo",
    help="The base url for SeqRepo REST API",
    show_default=True
)
def annotate_click(vcf_in: str, vcf_out: str, vrs_file: str,
                   seqrepo_dp_type: SeqRepoProxyType, seqrepo_root_dir: str,
                   seqrepo_base_url: str) -> None:
    """Annotate VCF file via click

    Example arguments:

    --vcf_in input.vcf.gz --vcf_out output.vcf.gz --vrs_file vrs_objects.pkl
    """
    annotator = VCFAnnotator(seqrepo_dp_type, seqrepo_base_url, seqrepo_root_dir)
    start = timer()
    annotator.annotate(vcf_in, vcf_out, vrs_file)
    end = timer()
    click.echo(f"VCF Annotator finished in {(end - start):.5f} seconds")

class VCFAnnotator:
    """Provides utility for annotating VCF's with VRS Allele IDs.
    VCF's are read using pysam and stored as pysam objects.
    Alleles are translated into VRS Allele IDs using VRS-Python Translator.
    """

    def __init__(self, seqrepo_dp_type: SeqRepoProxyType = SeqRepoProxyType.LOCAL,
                 seqrepo_base_url: str = "http://localhost:5000/seqrepo",
                 seqrepo_root_dir: str = "/usr/local/share/seqrepo/latest") -> None:
        """Initialize the VCFAnnotator class

        :param SeqRepoProxyType seqrepo_dp_type: The type of SeqRepo Data Proxy to use
        :param str seqrepo_base_url: The base url for SeqRepo REST API
        :param str seqrepo_root_dir: The root directory for the local SeqRepo instance
        """
        if seqrepo_dp_type == SeqRepoProxyType.LOCAL:
            self.dp = SeqRepoDataProxy(SeqRepo(seqrepo_root_dir))
        else:
            self.dp = SeqRepoRESTDataProxy(seqrepo_base_url)
        self.tlr = Translator(self.dp)

    def annotate(self, vcf_in: str, vcf_out: str, vrs_file: str) -> None:
        """Annotates an input VCF file with VRS Allele IDs & creates a pickle file
        containing the vrs object information.

        :param str vcf_in: The path for the input VCF file to annotate
        :param str vcf_out: The path for the output VCF file
        :param str vrs_file: The path for the output VCF pickle file
        """
        INFO_FIELD_ID = "VRS_Allele"  # pylint: disable=invalid-name
        vrs_data = {}
        vcf_in = pysam.VariantFile(filename=vcf_in)  # pylint: disable=no-member
        vcf_in.header.info.add(INFO_FIELD_ID, "1", "String", "vrs")
        vcf_out = pysam.VariantFile(vcf_out, "w", header=vcf_in.header)  # pylint: disable=no-member
        vrs_pickle_out = open(vrs_file, "wb")

        for record in vcf_in:
            vrs_allele_ids = self._record_digests(record, vrs_data)
            record.info[INFO_FIELD_ID] = ",".join(vrs_allele_ids)
            vcf_out.write(record)
        pickle.dump(vrs_data, vrs_pickle_out)


        vrs_pickle_out.close()
        vcf_in.close()
        vcf_out.close()

    def _get_vrs_object(self, vcf_coords: str, vrs_data: Dict,
                        vrs_allele_ids: List[str],
                        vrs_data_key: Optional[str] = None) -> None:
        """Get VRS Object given `vcf_coords`. `vrs_data` and `vrs_allele_ids` will
        be mutated.

        :param str vcf_coords: Allele to get VRS object for. Format is chr-pos-ref-alt
        :param Dict vrs_data: Dictionary containing the VRS object information for the
            VCF
        :param List[str] vrs_allele_ids: List containing the VRS Allele IDs
        :param Optional[str] vrs_data_key: The key to update in `vrs_data`. If not
            provided, will use `vcf_coords` as the key.
        """
        vrs_obj = self.tlr.translate_from(vcf_coords, "gnomad")
        key = vrs_data_key if vrs_data_key else vcf_coords
        vrs_data[key] = str(vrs_obj.as_dict())
        vrs_allele_ids.append(vrs_obj._id._value)

    def _record_digests(self, record: pysam.VariantRecord, vrs_data: Dict) -> List[str]:
        """Mutate `vrs_data` with VRS object information and returning a list of VRS
        Allele IDs

        :param pysam.VariantRecord record: A row in the VCF file
        :param Dict vrs_data: Dictionary containing the VRS object information for the
            VCF
        :return List[str] vrs_allele_ids: List containing the VRS Allele IDs
        """
        vrs_allele_ids = []

        # Get VRS data for reference allele
        gnomad_loc = f"{record.chrom}-{record.pos}"
        reference_allele = f"{gnomad_loc}-{record.ref}-{record.ref}"
        self._get_vrs_object(reference_allele, vrs_data, vrs_allele_ids)

        # Get VRS data for alts
        alts = record.alts or []
        alleles = [f"{gnomad_loc}-{record.ref}-{a}" for a in [*alts]]
        data = f"{record.chrom}\t{record.pos}\t{record.id}\t{record.ref}\t{record.alts}"
        for allele in alleles:
            if "*" in allele:
                vrs_allele_ids.append("")
            else:
                self._get_vrs_object(allele, vrs_data, vrs_allele_ids,
                                     vrs_data_key=data)
        return vrs_allele_ids


if __name__ == "__main__":
    # python3 -m src.ga4gh.vrs.extras.vcf_annotation --vcf_in input.vcf.gz \
    #    --vcf_out output.vcf.gz --vrs_file vrs_objects.pkl
    annotate_click()  # pylint: disable=no-value-for-parameter
