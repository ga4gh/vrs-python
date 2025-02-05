"""Annotate VCFs with VRS identifiers and attributes."""

import logging
import pickle
from enum import Enum
from pathlib import Path

import pysam
from biocommons.seqrepo import SeqRepo
from pydantic import ValidationError

from ga4gh.core.identifiers import (
    VrsObjectIdentifierIs,
    use_ga4gh_compute_identifier_when,
)
from ga4gh.vrs.dataproxy import (
    DataProxyValidationError,
    SeqRepoDataProxy,
    SeqRepoRESTDataProxy,
)
from ga4gh.vrs.extras.translator import AlleleTranslator

_logger = logging.getLogger(__name__)


class VCFAnnotatorError(Exception):
    """Custom exceptions for VCF Annotator tool"""


class SeqRepoProxyType(str, Enum):
    """Define constraints for SeqRepo Data Proxy types"""

    LOCAL = "local"
    REST = "rest"


class VCFAnnotator:
    """Annotate VCFs with VRS allele IDs.

    Uses pysam to read, store, and (optionally) output VCFs. Alleles are translated
    into VRS IDs using the VRS-Python translator class.
    """

    # Field names for VCF
    VRS_ALLELE_IDS_FIELD = "VRS_Allele_IDs"
    VRS_STARTS_FIELD = "VRS_Starts"
    VRS_ENDS_FIELD = "VRS_Ends"
    VRS_STATES_FIELD = "VRS_States"
    VRS_ERROR_FIELD = "VRS_Error"
    # VCF character escape map
    VCF_ESCAPE_MAP = [  # noqa: RUF012
        ("%", "%25"),
        (";", "%3B"),
        (",", "%2C"),
        ("\r", "%0D"),
        ("\n", "%0A"),
        ("\t", "%09"),
    ]

    def __init__(
        self,
        seqrepo_dp_type: SeqRepoProxyType = SeqRepoProxyType.LOCAL,
        seqrepo_base_url: str = "http://localhost:5000/seqrepo",
        seqrepo_root_dir: str = "/usr/local/share/seqrepo/latest",
    ) -> None:
        """Initialize the VCFAnnotator class.

        :param seqrepo_dp_type: The type of SeqRepo Data Proxy to use
            (i.e., local vs REST)
        :param seqrepo_base_url: The base url for SeqRepo REST API
        :param seqrepo_root_dir: The root directory for the local SeqRepo instance
        """
        if seqrepo_dp_type == SeqRepoProxyType.LOCAL:
            self.dp = SeqRepoDataProxy(SeqRepo(seqrepo_root_dir))
        else:
            self.dp = SeqRepoRESTDataProxy(seqrepo_base_url)
        self.tlr = AlleleTranslator(self.dp)

    @use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.MISSING)
    def annotate(
        self,
        input_vcf_path: Path,
        output_vcf_path: Path | None = None,
        output_pkl_path: Path | None = None,
        vrs_attributes: bool = False,
        assembly: str = "GRCh38",
        compute_for_ref: bool = True,
        require_validation: bool = True,
    ) -> None:
        """Given a VCF, produce an output VCF annotated with VRS allele IDs, and/or
        a pickle file containing the full VRS objects.

        :param input_vcf_path: Location of input VCF
        :param output_vcf_path: The path for the output VCF file
        :param output_pkl_path: The path for the output VCF pickle file
        :param vrs_attributes: If `True`, include VRS_Start, VRS_End, VRS_State
            properties in the VCF INFO field. If `False` will not include these
            properties. Only used if `vcf_out` is defined.
        :param assembly: The assembly used in `vcf_in` data
        :param compute_for_ref: If true, compute VRS IDs for the reference allele
        :param require_validation: If `True`, validation checks (i.e., REF value
            matches expected REF at given location) must pass in order to return a VRS
            object for a record. If `False` then VRS object will be returned even if
            validation checks fail, although all instances of failed validation are
            logged as warnings regardless.
        """
        if not any((output_vcf_path, output_pkl_path)):
            msg = "Must provide one of: `vcf_out` or `vrs_pickle_out`"
            raise VCFAnnotatorError(msg)

        info_field_num = "R" if compute_for_ref else "A"
        info_field_desc = "REF and ALT" if compute_for_ref else "ALT"

        vcf = pysam.VariantFile(filename=str(input_vcf_path.absolute()))
        vcf.header.info.add(
            self.VRS_ALLELE_IDS_FIELD,
            info_field_num,
            "String",
            (
                "The computed identifiers for the GA4GH VRS Alleles corresponding to the "
                f"GT indexes of the {info_field_desc} alleles"
            ),
        )
        vcf.header.info.add(
            self.VRS_ERROR_FIELD,
            ".",
            "String",
            ("If an error occurred computing a VRS Identifier, the error message"),
        )

        if vrs_attributes:
            vcf.header.info.add(
                self.VRS_STARTS_FIELD,
                info_field_num,
                "String",
                (
                    "Interresidue coordinates used as the location starts for the GA4GH "
                    f"VRS Alleles corresponding to the GT indexes of the {info_field_desc} alleles"
                ),
            )
            vcf.header.info.add(
                self.VRS_ENDS_FIELD,
                info_field_num,
                "String",
                (
                    "Interresidue coordinates used as the location ends for the GA4GH VRS "
                    f"Alleles corresponding to the GT indexes of the {info_field_desc} alleles"
                ),
            )
            vcf.header.info.add(
                self.VRS_STATES_FIELD,
                info_field_num,
                "String",
                (
                    "The literal sequence states used for the GA4GH VRS Alleles "
                    f"corresponding to the GT indexes of the {info_field_desc} alleles"
                ),
            )

        vcf_out = (
            pysam.VariantFile(str(output_vcf_path.absolute()), "w", header=vcf.header)
            if output_vcf_path
            else None
        )

        # only retain raw data if dumping to pkl
        vrs_data = {} if output_pkl_path else None
        for record in vcf:
            if vcf_out:
                additional_info_fields = [self.VRS_ALLELE_IDS_FIELD]
                if vrs_attributes:
                    additional_info_fields += [
                        self.VRS_STARTS_FIELD,
                        self.VRS_ENDS_FIELD,
                        self.VRS_STATES_FIELD,
                    ]
            else:
                # no INFO field names need to be designated if not producing an annotated VCF
                additional_info_fields = []
            try:
                vrs_field_data = self._get_vrs_data(
                    record,
                    vrs_data,
                    assembly,
                    additional_info_fields,
                    vrs_attributes=vrs_attributes,
                    compute_for_ref=compute_for_ref,
                    require_validation=require_validation,
                )
            except Exception as ex:
                _logger.exception("VRS error on %s-%s", record.chrom, record.pos)
                err_msg = f"{ex}" or f"{type(ex)}"
                for search_repl in VCFAnnotator.VCF_ESCAPE_MAP:
                    err_msg = err_msg.replace(search_repl[0], search_repl[1])
                additional_info_fields = [self.VRS_ERROR_FIELD]
                vrs_field_data = {self.VRS_ERROR_FIELD: [err_msg]}

            _logger.debug(
                "VCF record %s-%s generated vrs_field_data %s",
                record.chrom,
                record.pos,
                vrs_field_data,
            )

            if output_vcf_path and vcf_out:
                for k in additional_info_fields:
                    record.info[k] = [value or "." for value in vrs_field_data[k]]
                vcf_out.write(record)

        vcf.close()

        if vcf_out:
            vcf_out.close()

        if output_pkl_path:
            with output_pkl_path.open("wb") as wf:
                pickle.dump(vrs_data, wf)

    def _get_vrs_object(
        self,
        vcf_coords: str,
        vrs_data: dict | None,
        vrs_field_data: dict,
        assembly: str,
        vrs_data_key: str | None = None,
        vrs_attributes: bool = False,
        require_validation: bool = True,
    ) -> None:
        """Get VRS object given `vcf_coords`. `vrs_data` and `vrs_field_data` will
        be mutated.

        :param vcf_coords: Allele to get VRS object for. Format is chr-pos-ref-alt
        :param vrs_data: Dictionary containing the VRS object information for the VCF
        :param vrs_field_data: If `output_vcf`, keys are VRS Fields and values are list
            of VRS data. Empty otherwise.
        :param assembly: The assembly used in `vcf_coords`
        :param vrs_data_key: The key to update in `vrs_data`. If not provided, will use
            `vcf_coords` as the key.
        :param vrs_attributes: If `True` will include VRS_Start, VRS_End, VRS_State
            fields in the INFO field. If `False` will not include these fields. Only
            used if `output_vcf` set to `True`.
        :param require_validation: If `True` then validation checks must pass in order
            to return a VRS object. If `False` then VRS object will be returned even if
            validation checks fail. Defaults to `True`.
        """
        try:
            vrs_obj = self.tlr._from_gnomad(  # noqa: SLF001
                vcf_coords,
                assembly_name=assembly,
                require_validation=require_validation,
            )
        except (ValidationError, DataProxyValidationError):
            vrs_obj = None
            _logger.exception(
                "ValidationError when translating %s from gnomad", vcf_coords
            )
            raise
        except KeyError:
            vrs_obj = None
            _logger.exception("KeyError when translating %s from gnomad", vcf_coords)
            raise
        except AssertionError:
            vrs_obj = None
            _logger.exception(
                "AssertionError when translating %s from gnomad", vcf_coords
            )
            raise
        except Exception:
            vrs_obj = None
            _logger.exception(
                "Unhandled Exception when translating %s from gnomad", vcf_coords
            )
            raise
        else:
            if not vrs_obj:
                _logger.debug(
                    "None was returned when translating %s from gnomad", vcf_coords
                )

        if vrs_data and vrs_obj:
            key = vrs_data_key if vrs_data_key else vcf_coords
            vrs_data[key] = str(vrs_obj.model_dump(exclude_none=True))

        if vrs_field_data:
            allele_id = vrs_obj.id if vrs_obj else ""
            vrs_field_data[self.VRS_ALLELE_IDS_FIELD].append(allele_id)

            if vrs_attributes:
                if vrs_obj:
                    start = str(vrs_obj.location.start)
                    end = str(vrs_obj.location.end)
                    alt = (
                        str(vrs_obj.state.sequence.root)
                        if vrs_obj.state.sequence
                        else ""
                    )
                else:
                    start = end = alt = ""

                vrs_field_data[self.VRS_STARTS_FIELD].append(start)
                vrs_field_data[self.VRS_ENDS_FIELD].append(end)
                vrs_field_data[self.VRS_STATES_FIELD].append(alt)

    def _get_vrs_data(
        self,
        record: pysam.VariantRecord,
        vrs_data: dict | None,
        assembly: str,
        additional_info_fields: list[str],
        vrs_attributes: bool = False,
        compute_for_ref: bool = True,
        require_validation: bool = True,
    ) -> dict:
        """Get VRS data for record's reference and alt alleles.

        :param record: A row in the VCF file
        :param vrs_data: Dictionary containing the VRS object information for the VCF.
            Will be mutated if `output_pickle = True`
        :param assembly: The assembly used in `record`
        :param additional_info_fields: Additional VRS fields to add in INFO field
        :param vrs_attributes: If `True` will include VRS_Start, VRS_End, VRS_State
            fields in the INFO field. If `False` will not include these fields. Only
            used if `output_vcf` set to `True`.
        :param compute_for_ref: If true, compute VRS IDs for the reference allele
        :param require_validation: If `True` then validation checks must pass in
            order to return a VRS object. A `DataProxyValidationError` will be raised if
            validation checks fail. If `False` then VRS object will be returned even if
            validation checks fail. Defaults to `True`.
        :return: A dictionary mapping VRS-related INFO fields to lists of associated
            values. Will be empty if `create_annotated_vcf` is false.
        """
        vrs_field_data = {field: [] for field in additional_info_fields}

        # Get VRS data for reference allele
        gnomad_loc = f"{record.chrom}-{record.pos}"
        if compute_for_ref:
            reference_allele = f"{gnomad_loc}-{record.ref}-{record.ref}"
            self._get_vrs_object(
                reference_allele,
                vrs_data,
                vrs_field_data,
                assembly,
                vrs_attributes=vrs_attributes,
                require_validation=require_validation,
            )

        # Get VRS data for alts
        alts = record.alts or []
        alleles = [f"{gnomad_loc}-{record.ref}-{a}" for a in [*alts]]
        data = f"{record.chrom}\t{record.pos}\t{record.ref}\t{record.alts}"
        for allele in alleles:
            if "*" in allele:
                _logger.debug("Star allele found: %s", allele)
                for field in additional_info_fields:
                    vrs_field_data[field].append("")
            else:
                self._get_vrs_object(
                    allele,
                    vrs_data,
                    vrs_field_data,
                    assembly,
                    vrs_data_key=data,
                    vrs_attributes=vrs_attributes,
                    require_validation=require_validation,
                )

        return vrs_field_data
