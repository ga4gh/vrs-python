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
from ga4gh.vrs.extras.translator import Translator, ValidationError as TranslatorValidationError



_logger = logging.getLogger(__name__)
_logger.setLevel(logging.DEBUG)


class VCFAnnotatorException(Exception):
    """Custom exceptions for VCF Annotator tool"""


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
    "--vrs_attributes",
    is_flag=True,
    default=False,
    help="Will include VRS_Start, VRS_End, VRS_Alt fields in the INFO field.",
    show_default=True
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
    vrs_attributes: bool, seqrepo_dp_type: SeqRepoProxyType, seqrepo_root_dir: str,
    seqrepo_base_url: str, assembly: str
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
    annotator.annotate(vcf_in, vcf_out, vrs_pickle_out, vrs_attributes, assembly)
    end = timer()
    msg = f"VCF Annotator finished in {(end - start):.5f} seconds"
    _logger.info(msg)
    click.echo(msg)

class VCFAnnotator:  # pylint: disable=too-few-public-methods
    """Provides utility for annotating VCF's with VRS Allele IDs.
    VCF's are read using pysam and stored as pysam objects.
    Alleles are translated into VRS Allele IDs using VRS-Python Translator.
    """

    # Field names for VCF
    VRS_ALLELE_FIELD = "VRS_Allele"
    VRS_START_FIELD = "VRS_Start"
    VRS_END_FIELD = "VRS_End"
    VRS_ALT_FIELD = "VRS_Alt"

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

    def annotate(  # pylint: disable=too-many-arguments,too-many-locals
        self, vcf_in: str, vcf_out: Optional[str] = None,
        vrs_pickle_out: Optional[str] = None, vrs_attributes: bool = False,
        assembly: str = "GRCh38"
    ) -> None:
        """Annotates an input VCF file with VRS Allele IDs & creates a pickle file
        containing the vrs object information.

        :param str vcf_in: The path for the input VCF file to annotate
        :param Optional[str] vcf_out: The path for the output VCF file
        :param Optional[str] vrs_pickle_out: The path for the output VCF pickle file
        :param bool vrs_attributes: If `True` will include VRS_Start, VRS_End,
            VRS_Alt fields in the INFO field. If `False` will not include these fields.
            Only used if `vcf_out` is provided.
        :param str assembly: The assembly used in `vcf_in` data
        """
        if not any((vcf_out, vrs_pickle_out)):
            raise VCFAnnotatorException(
                "Must provide one of: `vcf_out` or `vrs_pickle_out`")

        vrs_data = {}
        vcf_in = pysam.VariantFile(filename=vcf_in)  # pylint: disable=no-member
        vcf_in.header.info.add(
            self.VRS_ALLELE_FIELD, 1, "String",
            "VRS Allele IDs, separated by commas, associated to REF and ALT fields."
        )

        additional_info_fields = [self.VRS_ALLELE_FIELD]
        if vrs_attributes:
            vcf_in.header.info.add(
                self.VRS_START_FIELD, 1, "String",
                "VRS Allele start positions (inter-residue) associated to REF and ALT fields."
            )
            vcf_in.header.info.add(
                self.VRS_END_FIELD, 1, "String",
                "VRS Allele end positions (inter-residue) associated to REF and ALT fields."
            )
            vcf_in.header.info.add(
                self.VRS_ALT_FIELD, 1, "String",
                "VRS Allele alt sequences associated to REF and ALT fields."
            )
            additional_info_fields += [self.VRS_START_FIELD, self.VRS_END_FIELD, self.VRS_ALT_FIELD]

        if vcf_out:
            vcf_out = pysam.VariantFile(vcf_out, "w", header=vcf_in.header)  # pylint: disable=no-member

        output_vcf = bool(vcf_out)
        output_pickle = bool(vrs_pickle_out)

        for record in vcf_in:
            vrs_field_data = self._get_vrs_data(
                record, vrs_data, assembly, additional_info_fields, vrs_attributes,
                output_pickle, output_vcf
            )

            if output_vcf:
                for k in additional_info_fields:
                    record.info[k] = ",".join(vrs_field_data[k])
                vcf_out.write(record)

        vcf_in.close()

        if output_vcf:
            vcf_out.close()

        if vrs_pickle_out:
            with open(vrs_pickle_out, "wb") as wf:
                pickle.dump(vrs_data, wf)

    def _get_vrs_object(  # pylint: disable=too-many-arguments,too-many-locals
        self, vcf_coords: str, vrs_data: Dict, vrs_field_data: Dict, assembly: str,
        vrs_data_key: Optional[str] = None, output_pickle: bool = True,
        output_vcf: bool = False, vrs_attributes: bool = False
    ) -> None:
        """Get VRS Object given `vcf_coords`. `vrs_data` and `vrs_field_data` will
        be mutated.

        :param str vcf_coords: Allele to get VRS object for. Format is chr-pos-ref-alt
        :param Dict vrs_data: Dictionary containing the VRS object information for the
            VCF
        :param Dict vrs_field_data: If `output_vcf`, will keys will be VRS Fields and
            values will be list of VRS data. Else, will be an empty dictionary
        :param str assembly: The assembly used in `vcf_coords`
        :param Optional[str] vrs_data_key: The key to update in `vrs_data`. If not
            provided, will use `vcf_coords` as the key.
        :param bool output_pickle: `True` if VRS pickle file will be output.
            `False` otherwise.
        :param bool output_vcf: `True` if annotated VCF file will be output.
            `False` otherwise.
        :param bool vrs_attributes: If `True` will include VRS_Start, VRS_End,
            VRS_Alt fields in the INFO field. If `False` will not include these fields.
            Only used if `output_vcf` set to `True`.
        """
        try:
            vrs_obj = self.tlr._from_gnomad(vcf_coords, assembly_name=assembly)
        except (ValidationError, TranslatorValidationError) as e:
            vrs_obj = None
            _logger.error("ValidationError when translating %s from gnomad: %s", vcf_coords, str(e))
        except KeyError as e:
            vrs_obj = None
            _logger.error("KeyError when translating %s from gnomad: %s", vcf_coords, str(e))
        except AssertionError as e:
            vrs_obj = None
            _logger.error("AssertionError when translating %s from gnomad: %s", vcf_coords, str(e))
        except Exception as e:  # pylint: disable=broad-except
            vrs_obj = None
            _logger.error("Unhandled Exception when translating %s from gnomad: %s", vcf_coords, str(e))
        else:
            if not vrs_obj:
                _logger.debug("None was returned when translating %s from gnomad", vcf_coords)

        if output_pickle and vrs_obj:
            key = vrs_data_key if vrs_data_key else vcf_coords
            vrs_data[key] = str(vrs_obj.as_dict())

        if output_vcf:
            allele_id = vrs_obj._id._value if vrs_obj else ""
            vrs_field_data[self.VRS_ALLELE_FIELD].append(allele_id)

            if vrs_attributes:
                start = str(vrs_obj.location.interval.start.value) if vrs_obj else ""
                end = str(vrs_obj.location.interval.end.value) if vrs_obj else ""
                alt = str(vrs_obj.state.sequence) if vrs_obj else ""

                vrs_field_data[self.VRS_START_FIELD].append(start)
                vrs_field_data[self.VRS_END_FIELD].append(end)
                vrs_field_data[self.VRS_ALT_FIELD].append(alt)

    def _get_vrs_data(  # pylint: disable=too-many-arguments,too-many-locals
        self, record: pysam.VariantRecord, vrs_data: Dict, assembly: str,  # pylint: disable=no-member
        additional_info_fields: List[str], vrs_attributes: bool = False,
        output_pickle: bool = True, output_vcf: bool = True
    ) -> Dict:
        """Get VRS data for record's reference and alt alleles.

        :param pysam.VariantRecord record: A row in the VCF file
        :param Dict vrs_data: Dictionary containing the VRS object information for the
            VCF. Will be mutated if `output_pickle = True`
        :param str assembly: The assembly used in `record`
        :param List[str] additional_info_fields: Additional VRS fields to add in INFO
            field
        :param bool vrs_attributes: If `True` will include VRS_Start, VRS_End,
            VRS_Alt fields in the INFO field. If `False` will not include these fields.
            Only used if `output_vcf` set to `True`.
        :param bool output_pickle: `True` if VRS pickle file will be output.
            `False` otherwise.
        :param bool output_vcf: `True` if annotated VCF file will be output.
            `False` otherwise.
        :return: If `output_vcf = True`, a dictionary containing VRS Fields and list
            of associated values. If `output_vcf = False`, an empty dictionary will be
            returned.
        """
        vrs_field_data = {field: [] for field in additional_info_fields} if output_vcf else {}

        # Get VRS data for reference allele
        gnomad_loc = f"{record.chrom}-{record.pos}"
        reference_allele = f"{gnomad_loc}-{record.ref}-{record.ref}"
        self._get_vrs_object(
            reference_allele, vrs_data, vrs_field_data, assembly,
            output_pickle=output_pickle, output_vcf=output_vcf,
            vrs_attributes=vrs_attributes
        )

        # Get VRS data for alts
        alts = record.alts or []
        alleles = [f"{gnomad_loc}-{record.ref}-{a}" for a in [*alts]]
        data = f"{record.chrom}\t{record.pos}\t{record.ref}\t{record.alts}"
        for allele in alleles:
            if "*" in allele:
                _logger.debug("Star allele found: %s", allele)
                if output_vcf:
                    for field in additional_info_fields:
                        vrs_field_data[field].append("")
            else:
                self._get_vrs_object(
                    allele, vrs_data, vrs_field_data, assembly, vrs_data_key=data,
                    output_pickle=output_pickle, output_vcf=output_vcf,
                    vrs_attributes=vrs_attributes
                )

        return vrs_field_data


if __name__ == "__main__":
    # python3 -m src.ga4gh.vrs.extras.vcf_annotation --vcf_in input.vcf.gz \
    #    --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl
    annotate_click()  # pylint: disable=no-value-for-parameter
