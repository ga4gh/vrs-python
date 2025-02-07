"""Annotate VCFs with VRS identifiers and attributes."""

import abc
import logging
from enum import Enum
from pathlib import Path

import pysam

from ga4gh.core.identifiers import (
    VrsObjectIdentifierIs,
    use_ga4gh_compute_identifier_when,
)
from ga4gh.vrs.dataproxy import _DataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.models import Allele

_logger = logging.getLogger(__name__)


class VCFAnnotatorError(Exception):
    """Custom exceptions for VCF Annotator tool"""


class FieldName(str, Enum):
    """Define VCF field names for VRS annotations"""

    IDS_FIELD = "VRS_Allele_IDs"
    STARTS_FIELD = "VRS_Starts"
    ENDS_FIELD = "VRS_Ends"
    STATES_FIELD = "VRS_States"
    ERROR_FIELD = "VRS_Error"


# VCF character escape map
VCF_ESCAPE_MAP = str.maketrans(
    {
        "%": "%25",
        ";": "%3B",
        ",": "%2C",
        "\r": "%0D",
        "\n": "%0A",
    }
)


class AbstractVcfAnnotator(abc.ABC):
    """Abstract class for VCF annotation with VRS."""

    _collect_alleles: bool = False

    def __init__(self, data_proxy: _DataProxy) -> None:
        """Initialize the VCFAnnotator class.

        :param data_proxy: GA4GH sequence dataproxy instance.
        """
        self.data_proxy = data_proxy
        self.tlr = AlleleTranslator(self.data_proxy)

    @abc.abstractmethod
    def raise_for_output_args(self, output_vcf_path: Path | None, **kwargs) -> None:
        """Raise an exception if no output appears to be configured or declared.

        Child classes should implement this to ensure that any other outputs are checked
        on top of an annotated VCF output.

        :raise VCFAnnotatorError: if no output args are shown
        """

    def _update_vcf_header(
        self, vcf: pysam.VariantFile, incl_ref_allele: bool, incl_vrs_attrs: bool
    ) -> None:
        """Add new fields to VCF header

        :param vcf: pysam VCF object to annotate
        :param incl_ref_allele: whether VRS alleles will be calculated for REFs
        :param incl_vrs_attrs: whether INFO properties should be defined for VRS attributes
            (normalized coordinates/state)
        """
        info_field_num = "R" if incl_ref_allele else "A"
        info_field_desc = "REF and ALT" if incl_ref_allele else "ALT"

        vcf.header.info.add(
            FieldName.IDS_FIELD.value,
            info_field_num,
            "String",
            (
                "The computed identifiers for the GA4GH VRS Alleles corresponding to the "
                f"GT indexes of the {info_field_desc} alleles"
            ),
        )
        vcf.header.info.add(
            FieldName.ERROR_FIELD.value,
            ".",
            "String",
            ("If an error occurred computing a VRS Identifier, the error message"),
        )

        if incl_vrs_attrs:
            vcf.header.info.add(
                FieldName.STARTS_FIELD.value,
                info_field_num,
                "String",
                (
                    "Interresidue coordinates used as the location starts for the GA4GH "
                    f"VRS Alleles corresponding to the GT indexes of the {info_field_desc} alleles"
                ),
            )
            vcf.header.info.add(
                FieldName.ENDS_FIELD.value,
                info_field_num,
                "String",
                (
                    "Interresidue coordinates used as the location ends for the GA4GH VRS "
                    f"Alleles corresponding to the GT indexes of the {info_field_desc} alleles"
                ),
            )
            vcf.header.info.add(
                FieldName.STATES_FIELD.value,
                info_field_num,
                "String",
                (
                    "The literal sequence states used for the GA4GH VRS Alleles "
                    f"corresponding to the GT indexes of the {info_field_desc} alleles"
                ),
            )

    @use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.MISSING)
    def annotate(
        self,
        input_vcf_path: Path,
        output_vcf_path: Path | None = None,
        vrs_attributes: bool = False,
        assembly: str = "GRCh38",
        compute_for_ref: bool = True,
        require_validation: bool = True,
        **kwargs,
    ) -> None:
        """Given a VCF, produce an output VCF annotated with VRS allele IDs, and/or
        a pickle file containing the full VRS objects.

        :param input_vcf_path: Location of input VCF
        :param output_vcf_path: The path for the output VCF file
        :param output_pkl_path: The path for the output VCF pickle file
        :param vrs_attributes: If `True`, include VRS_Start, VRS_End, VRS_State
            properties in the VCF INFO field. If `False` will not include these
            properties. Only used if `output_vcf_path` is defined.
        :param assembly: The assembly used in `input_vcf_path` data
        :param compute_for_ref: If true, compute VRS IDs for the reference allele
        :param require_validation: If `True`, validation checks (i.e., REF value
            matches expected REF at given location) must pass in order to return a VRS
            object for a record. If `False` then VRS object will be returned even if
            validation checks fail, although all instances of failed validation are
            logged as warnings regardless.
        :raise VCFAnnotatorError: if no output formats are selected
        """
        self.raise_for_output_args(output_vcf_path, **kwargs)

        vcf = pysam.VariantFile(filename=str(input_vcf_path.absolute()))
        if output_vcf_path:
            self._update_vcf_header(vcf, compute_for_ref, vrs_attributes)
            vcf_out = pysam.VariantFile(
                str(output_vcf_path.absolute()), "w", header=vcf.header
            )
        else:
            vcf_out = None

        allele_collection = [] if self._collect_alleles else None
        for record in vcf:
            if vcf_out:
                additional_info_fields = [FieldName.IDS_FIELD]
                if vrs_attributes:
                    additional_info_fields += [
                        FieldName.STARTS_FIELD,
                        FieldName.ENDS_FIELD,
                        FieldName.STATES_FIELD,
                    ]
            else:
                # no INFO field names need to be designated if not producing an annotated VCF
                additional_info_fields = []
            try:
                vrs_field_data = self._get_vrs_data(
                    record,
                    allele_collection,
                    assembly,
                    additional_info_fields,
                    vrs_attributes=vrs_attributes,
                    compute_for_ref=compute_for_ref,
                    require_validation=require_validation,
                    **kwargs,
                )
            except Exception as ex:
                _logger.exception("VRS error on %s-%s", record.chrom, record.pos)
                err_msg = f"{ex}" or f"{type(ex)}"
                err_msg = err_msg.translate(VCF_ESCAPE_MAP)
                additional_info_fields = [FieldName.ERROR_FIELD]
                vrs_field_data = {FieldName.ERROR_FIELD.value: [err_msg]}

            _logger.debug(
                "VCF record %s-%s generated vrs_field_data %s",
                record.chrom,
                record.pos,
                vrs_field_data,
            )

            if output_vcf_path and vcf_out:
                for k in additional_info_fields:
                    record.info[k.value] = [
                        value or "." for value in vrs_field_data[k.value]
                    ]
                vcf_out.write(record)

        vcf.close()

        if vcf_out:
            vcf_out.close()

        self.on_vrs_object_collection(allele_collection, **kwargs)

    @abc.abstractmethod
    def on_vrs_object(self, vcf_coords: str, vrs_allele: Allele, **kwargs) -> Allele:
        """Perform side-effects (eg additional annotation or storage) on VRS alleles
        as they are constructed during VCF annotation.

        Reimplement in a child class to add custom logic. Otherwise, this method simply
        passes through ``vrs_allele`` without altering it further or storing it.

        :param vcf_coords: CHR-POS-REF-ALT from VCF for this allele
        :param vrs_allele: allele translated from coords
        :return: final VRS allele
        """

    @abc.abstractmethod
    def on_vrs_object_collection(
        self, vrs_alleles_collection: list[Allele] | None, **kwargs
    ) -> None:
        """Perform clean-up operations (eg file writing) on VRS objects collected
        during VCF ingestion.

        Reimplement in a child class to add custom logic. Otherwise, this method does
        nothing.

        :param vrs_alleles: VRS alleles constructed from ingested VCF
        """

    def _get_vrs_object(
        self,
        vcf_coords: str,
        allele_collection: list[Allele] | None,
        vrs_field_data: dict,
        assembly: str,
        vrs_attributes: bool = False,
        require_validation: bool = True,
        **kwargs,
    ) -> None:
        """Get VRS object given `vcf_coords`. `vrs_data` and `vrs_field_data` will
        be mutated.

        :param vcf_coords: Allele to get VRS object for. Format is chr-pos-ref-alt
        :param vrs_data: All constructed VRS objects. Can be `None` if no data dumps
            will be created.
        :param vrs_field_data: If `vrs_data`, keys are VRS Fields and values are list
            of VRS data. Empty dict otherwise.
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
        except Exception:
            vrs_obj = None
            _logger.exception(
                "Exception encountered during translation of variation: %s", vcf_coords
            )
            raise
        if vrs_obj is None:
            _logger.debug(
                "None was returned when translating %s from gnomad", vcf_coords
            )
        else:
            vrs_obj = self.on_vrs_object(vcf_coords, vrs_obj, **kwargs)

        if allele_collection and vrs_obj:
            allele_collection.append(vrs_obj)

        if vrs_field_data:
            allele_id = vrs_obj.id if vrs_obj else ""
            vrs_field_data[FieldName.IDS_FIELD].append(allele_id)

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

                vrs_field_data[FieldName.STARTS_FIELD].append(start)
                vrs_field_data[FieldName.ENDS_FIELD].append(end)
                vrs_field_data[FieldName.STATES_FIELD].append(alt)

    def _get_vrs_data(
        self,
        record: pysam.VariantRecord,
        allele_collection: list | None,
        assembly: str,
        additional_info_fields: list[FieldName],
        vrs_attributes: bool = False,
        compute_for_ref: bool = True,
        require_validation: bool = True,
        **kwargs,
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
        vrs_field_data = {field.value: [] for field in additional_info_fields}

        # Get VRS data for reference allele
        gnomad_loc = f"{record.chrom}-{record.pos}"
        if compute_for_ref:
            reference_allele = f"{gnomad_loc}-{record.ref}-{record.ref}"
            self._get_vrs_object(
                reference_allele,
                allele_collection,
                vrs_field_data,
                assembly,
                vrs_attributes=vrs_attributes,
                require_validation=require_validation,
                **kwargs,
            )

        # Get VRS data for alts
        alts = record.alts or []
        alleles = [f"{gnomad_loc}-{record.ref}-{a}" for a in [*alts]]
        data = f"{record.chrom}\t{record.pos}\t{record.ref}\t{record.alts}"
        for allele in alleles:
            if "*" in allele:
                _logger.debug("Star allele found: %s", allele)
                for field in additional_info_fields:
                    vrs_field_data[field.value].append("")
            else:
                self._get_vrs_object(
                    allele,
                    allele_collection,
                    vrs_field_data,
                    assembly,
                    vrs_data_key=data,
                    vrs_attributes=vrs_attributes,
                    require_validation=require_validation,
                )

        return vrs_field_data


class VcfAnnotator(AbstractVcfAnnotator):
    """Annotate VCFs with VRS allele IDs.

    Uses pysam to ingest a VCF and produce a copy annotated with VRS IDs
    (and, optionally, normalized allele attributes). Alleles are translated
    using the VRS-Python translator class.
    """

    def raise_for_output_args(self, output_vcf_path: Path | None, **kwargs) -> None:  # noqa: D102 ARG002
        if output_vcf_path is None:
            msg = "No VCF output location provided."
            raise VCFAnnotatorError(msg)

    def on_vrs_object(self, vcf_coords: str, vrs_allele: Allele, **kwargs) -> Allele:  # noqa: ARG002
        """Perform side-effects (eg additional annotation or storage) on VRS alleles
        as they are constructed during VCF annotation.

        For this basic implementation, no additional operations are performed.

        :param vcf_coords: CHR-POS-REF-ALT from VCF for this allele
        :param vrs_allele: allele translated from coords
        :return: final VRS allele
        """
        return vrs_allele

    def on_vrs_object_collection(
        self, vrs_alleles_collection: list[Allele] | None, **kwargs
    ) -> None:
        """Perform clean-up operations (eg file writing) on VRS objects collected
        during VCF ingestion.

        For this basic implementation, does nothing.

        :param vrs_alleles: VRS alleles constructed from ingested VCF
        """
