"""Module containing tools for annotating VCFs with VRS

Example of how to run from root of vrs-python directory:
python3 -m src.ga4gh.vrs.extras.vcf_annotation --vcf_in input.vcf.gz \
    --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl
"""
import logging
import pickle
from enum import Enum
from typing import Dict, List, Optional
from timeit import default_timer as timer

import click
import pysam
from biocommons.seqrepo import SeqRepo
from python_jsonschema_objects import ValidationError

from ga4gh.vrs.dataproxy import SeqRepoDataProxy, SeqRepoRESTDataProxy
from ga4gh.vrs.extras.translator import Translator



_logger = logging.getLogger(__name__)
_logger.setLevel(logging.DEBUG)


class VCFAnnotatorException(Exception):
    """Custom exceptions for VCF Annotator tool"""

    pass  # pylint: disable=unnecessary-pass


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
    required=False,
    type=str,
    help=("The path for the output VCF file. If not provided, must provide "
          "--vrs_pickle_out.")
)
@click.option(
    "--vrs_pickle_out",
    required=False,
    type=str,
    help=("The path for the output VCF pickle file. If not provided, must provide "
          "--vcf_out")
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
@click.option(
    "--assembly",
    required=False,
    default="GRCh38",
    show_default=True,
    help="The assembly that the `vcf_in` data uses.",
    type=str
)
def annotate_click(  # pylint: disable=too-many-arguments
    vcf_in: str, vcf_out: Optional[str], vrs_pickle_out: Optional[str],
    seqrepo_dp_type: SeqRepoProxyType, seqrepo_root_dir: str, seqrepo_base_url: str,
    assembly: str
) -> None:
    """Annotate VCF file via click

    Example arguments:

    --vcf_in input.vcf.gz --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl
    """
    annotator = VCFAnnotator(seqrepo_dp_type, seqrepo_base_url, seqrepo_root_dir)
    start = timer()
    msg = f"Annotating {vcf_in} with the VCF Annotator..."
    _logger.info(msg)
    click.echo(msg)
    annotator.annotate(vcf_in, vcf_out, vrs_pickle_out, assembly)
    end = timer()
    msg = f"VCF Annotator finished in {(end - start):.5f} seconds"
    _logger.info(msg)
    click.echo(msg)

class VCFAnnotator:  # pylint: disable=too-few-public-methods
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

    def annotate(
        self, vcf_in: str, vcf_out: Optional[str] = None,
        vrs_pickle_out: Optional[str] = None, assembly: str = "GRCh38"
    ) -> None:
        """Annotates an input VCF file with VRS Allele IDs & creates a pickle file
        containing the vrs object information.

        :param str vcf_in: The path for the input VCF file to annotate
        :param Optional[str] vcf_out: The path for the output VCF file
        :param Optional[str] vrs_pickle_out: The path for the output VCF pickle file
        :param str assembly: The assembly used in `vcf_in` data
        """
        if not any((vcf_out, vrs_pickle_out)):
            raise VCFAnnotatorException(
                "Must provide one of: `vcf_out` or `vrs_pickle_out`")

        INFO_FIELD_ID = "VRS_Allele"  # pylint: disable=invalid-name
        vrs_data = {}
        vcf_in = pysam.VariantFile(filename=vcf_in)  # pylint: disable=no-member
        vcf_in.header.info.add(INFO_FIELD_ID, "1", "String", "vrs")
        if vcf_out:
            vcf_out = pysam.VariantFile(vcf_out, "w", header=vcf_in.header)  # pylint: disable=no-member

        output_vcf = bool(vcf_out)
        output_pickle = bool(vrs_pickle_out)

        for record in vcf_in:
            vrs_allele_ids = self._record_digests(record, vrs_data, assembly,
                                                  output_pickle, output_vcf)
            if output_vcf:
                record.info[INFO_FIELD_ID] = ",".join(vrs_allele_ids)
                vcf_out.write(record)
        vcf_in.close()

        if output_vcf:
            vcf_out.close()

        if vrs_pickle_out:
            with open(vrs_pickle_out, "wb") as wf:
                pickle.dump(vrs_data, wf)

    def _get_vrs_object(  # pylint: disable=too-many-arguments
        self, vcf_coords: str, vrs_data: Dict, vrs_allele_ids: List[str], assembly: str,
        vrs_data_key: Optional[str] = None, output_pickle: bool = True,
        output_vcf: bool = False
    ) -> None:
        """Get VRS Object given `vcf_coords`. `vrs_data` and `vrs_allele_ids` will
        be mutated.

        :param str vcf_coords: Allele to get VRS object for. Format is chr-pos-ref-alt
        :param Dict vrs_data: Dictionary containing the VRS object information for the
            VCF
        :param List[str] vrs_allele_ids: List containing the VRS Allele IDs
        :param str assembly: The assembly used in `vcf_coords`
        :param Optional[str] vrs_data_key: The key to update in `vrs_data`. If not
            provided, will use `vcf_coords` as the key.
        :param bool output_pickle: `True` if VRS pickle file will be output.
            `False` otherwise.
        :param bool output_vcf: `True` if annotated VCF file will be output.
            `False` otherwise.
        """
        try:
            vrs_obj = self.tlr._from_gnomad(vcf_coords, assembly_name=assembly)  # pylint: disable=protected-access
        except ValidationError as e:
            _logger.error("ValidationError when translating %s from gnomad: %s", vcf_coords, str(e))
        except KeyError as e:
            _logger.error("KeyError when translating %s from gnomad: %s", vcf_coords, str(e))
        except AssertionError as e:
            _logger.error("AssertionError when translating %s from gnomad: %s", vcf_coords, str(e))
        except Exception as e:  # pylint: disable=broad-except
            _logger.error("Unhandled Exception when translating %s from gnomad: %s", vcf_coords, str(e))
        else:
            if vrs_obj:
                if output_pickle:
                    key = vrs_data_key if vrs_data_key else vcf_coords
                    vrs_data[key] = str(vrs_obj.as_dict())

                if output_vcf:
                    vrs_allele_ids.append(vrs_obj._id._value)  # pylint: disable=protected-access
            else:
                _logger.debug("None was returned when translating %s from gnomad", vcf_coords)

    def _record_digests(  # pylint: disable=too-many-arguments
        self, record: pysam.VariantRecord, vrs_data: Dict, assembly: str,  # pylint: disable=no-member
        output_pickle: bool = True, output_vcf: bool = True
    ) -> List[str]:
        """Get VRS data for record's reference and alt alleles.

        :param pysam.VariantRecord record: A row in the VCF file
        :param Dict vrs_data: Dictionary containing the VRS object information for the
            VCF. Will be mutated if `output_pickle = True`
        :param str assembly: The assembly used in `record`
        :param bool output_pickle: `True` if VRS pickle file will be output.
            `False` otherwise.
        :param bool output_vcf: `True` if annotated VCF file will be output.
            `False` otherwise.
        :return List[str] vrs_allele_ids: List containing the VRS Allele IDs.
            If `output_vcf = False`, an empty list will be returned.
        """
        vrs_allele_ids = []

        # Get VRS data for reference allele
        gnomad_loc = f"{record.chrom}-{record.pos}"
        reference_allele = f"{gnomad_loc}-{record.ref}-{record.ref}"
        self._get_vrs_object(reference_allele, vrs_data, vrs_allele_ids, assembly,
                             output_pickle=output_pickle,
                             output_vcf=output_vcf)

        # Get VRS data for alts
        alts = record.alts or []
        alleles = [f"{gnomad_loc}-{record.ref}-{a}" for a in [*alts]]
        data = f"{record.chrom}\t{record.pos}\t{record.ref}\t{record.alts}"
        for allele in alleles:
            if "*" in allele:
                _logger.debug("Star allele found: %s", allele)
                vrs_allele_ids.append("")
            else:
                self._get_vrs_object(allele, vrs_data, vrs_allele_ids, assembly,
                                     vrs_data_key=data, output_pickle=output_pickle,
                                     output_vcf=output_vcf)
        return vrs_allele_ids


if __name__ == "__main__":
    # python3 -m src.ga4gh.vrs.extras.vcf_annotation --vcf_in input.vcf.gz \
    #    --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl
    annotate_click()  # pylint: disable=no-value-for-parameter
