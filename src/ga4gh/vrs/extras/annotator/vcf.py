"""Annotate VCFs with VRS identifiers and attributes, and read VRS Alleles back from
annotated VCFs.
"""

import abc
import logging
import pickle
from collections.abc import Iterator
from enum import Enum
from pathlib import Path
from typing import Literal

import pysam
from pydantic import ValidationError

from ga4gh.core import ga4gh_identify
from ga4gh.core.identifiers import (
    VrsObjectIdentifierIs,
    use_ga4gh_compute_identifier_when,
)
from ga4gh.vrs import VRS_VERSION, VrsType, __version__
from ga4gh.vrs.dataproxy import _DataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.models import (
    Allele,
    LengthExpression,
    LiteralSequenceExpression,
    Range,
    ReferenceLengthExpression,
    SequenceLocation,
    SequenceReference,
)

_logger = logging.getLogger(__name__)


class VcfAnnotatorError(Exception):
    """Custom exceptions for VCF Annotator tool"""


class VcfAnnotatorArgsError(VcfAnnotatorError):
    """Raise for improper args passed to VCF annotator"""


class VcfAnnotatorAbortError(VcfAnnotatorError):
    """Raise to abort annotation of the entire VCF.

    By default, an exception raised while processing a VCF record is logged and
    recorded in that record's ``VRS_Error`` INFO field, and annotation continues
    with the next record. Child classes of ``AbstractVcfAnnotator`` can raise this
    exception (e.g. from ``on_vrs_object``) to signal an unrecoverable error, such as
    a lost database connection, that should instead stop annotation and propagate
    to the caller of ``annotate()``.
    """


class VcfReaderError(Exception):
    """Raise for VCFs that VRS Alleles can't be read back from"""


class FieldName(str, Enum):
    """Define VCF field names for VRS annotations"""

    IDS_FIELD = "VRS_Allele_IDs"
    STARTS_FIELD = "VRS_Starts"
    ENDS_FIELD = "VRS_Ends"
    STATES_FIELD = "VRS_States"
    LENGTHS_FIELD = "VRS_Lengths"
    REPEAT_SUBUNIT_LENGTHS_FIELD = "VRS_RepeatSubunitLengths"
    ERROR_FIELD = "VRS_Error"


# String-typed INFO fields where pysam incorrectly converts None → b""
# (empty bytes) instead of the VCF missing value ".". Integer/Float fields
# are unaffected because pysam uses proper BCF missing sentinels for those.
_STRING_FIELDS = frozenset(
    {
        FieldName.IDS_FIELD,
        FieldName.STATES_FIELD,
        FieldName.ERROR_FIELD,
    }
)


# INFO fields, in addition to ``FieldName.IDS_FIELD``, that ``VcfAnnotator`` adds when
# run with VRS attributes. ``VcfReader`` needs all of them to reconstruct VRS Alleles.
_VRS_ATTRIBUTE_FIELDS = (
    FieldName.STARTS_FIELD,
    FieldName.ENDS_FIELD,
    FieldName.STATES_FIELD,
    FieldName.LENGTHS_FIELD,
    FieldName.REPEAT_SUBUNIT_LENGTHS_FIELD,
)


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

# ReferenceLengthExpression .sequence values will be included in output VCF if
# length <= this value. This field is optional for RLE since it can be derived
# from the reference sequence. Set to None to always include the sequence.
RLE_SEQ_LIMIT = 50


def dump_alleles_to_pkl(alleles: list[Allele], output_pkl_path: Path) -> None:
    """Create pkl file of dictionary mapping VRS IDs to ingested alleles.

    :param alleles: all alleles constructed + annotated during VCF ingestion
    :param output_pkl_path: location to save PKL file to
    """
    allele_dict = {
        allele.id: allele.model_dump(exclude_none=True) for allele in alleles
    }
    with output_pkl_path.open("wb") as f:
        pickle.dump(allele_dict, f)


def dump_alleles_to_ndjson(
    allele_collection: list[Allele], output_ndjson_path: Path
) -> None:
    """Create NDJSON dump of all alleles ingested from VCF

    :param alleles: all alleles + constructed + annotated during VCF ingestion
    :param output_ndjson_path: location to save NDJSON file to
    """
    with output_ndjson_path.open("w") as f:
        for allele in allele_collection:
            f.write(allele.model_dump_json(exclude_none=True) + "\n")


class AbstractVcfAnnotator(abc.ABC):
    """Abstract class for VCF annotation with VRS.

    Child classes must implement all abstract methods, though some may be less relevant
    and can just be stubbed out (for example, if the desired side effect for each VRS
    object is something like a REST request, there may be no need for extra cleanup in
    ``on_vrs_object_collection``).

    Critically, implementations that intend to do something with all VRS objects at once
    following ingestion (e.g. writing a PKL dump) need to set the class variable
    ``collect_alleles`` to ``True``.
    """

    collect_alleles: bool = False

    def __init__(self, data_proxy: _DataProxy, **kwargs) -> None:  # noqa: ARG002
        """Initialize the VCFAnnotator class.

        :param data_proxy: GA4GH sequence dataproxy instance.
        """
        self.data_proxy = data_proxy
        self.tlr = AlleleTranslator(self.data_proxy)

    def should_collect_alleles(self, **kwargs) -> bool:  # noqa: ARG002
        """Determine whether allele aggregation is necessary.

        This method is called to initialize (or not) an initial allele collection for
        downstream use. By default, it returns the corresponding class variable, but
        implementing classes can choose to amend this based on provided annotation
        arguments.

        :return: ``True`` if alleles should be collected
        """
        return self.collect_alleles

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
                f"The computed identifiers for the GA4GH VRS Alleles corresponding to the GT indexes of the {info_field_desc} alleles [VRS version={VRS_VERSION};VRS-Python version={__version__}]"
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
                "Integer",
                (
                    "Interresidue coordinates used as the location starts for the GA4GH "
                    f"VRS Alleles corresponding to the GT indexes of the {info_field_desc} alleles"
                ),
            )
            vcf.header.info.add(
                FieldName.ENDS_FIELD.value,
                info_field_num,
                "Integer",
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
            vcf.header.info.add(
                FieldName.LENGTHS_FIELD.value,
                info_field_num,
                "Integer",
                (
                    "The length values from ReferenceLengthExpression states for the GA4GH VRS "
                    f"Alleles corresponding to the GT indexes of the {info_field_desc} alleles"
                ),
            )
            vcf.header.info.add(
                FieldName.REPEAT_SUBUNIT_LENGTHS_FIELD.value,
                info_field_num,
                "Integer",
                (
                    "The repeatSubunitLength values from ReferenceLengthExpression states for the GA4GH VRS "
                    f"Alleles corresponding to the GT indexes of the {info_field_desc} alleles"
                ),
            )

    @use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.MISSING)
    def annotate(
        self,
        input_vcf_path: Path | Literal["-"],
        output_vcf_path: Path | Literal["-"] | None = None,
        vrs_attributes: bool = False,
        assembly: str = "GRCh38",
        compute_for_ref: bool = True,
        require_validation: bool = True,
        **kwargs,
    ) -> None:
        """Given a VCF, produce an output VCF annotated with VRS allele IDs, and/or
        additional storage outputs as implemented in a specific child class.

        :param input_vcf_path: Location of input VCF (or "-" to read from stdin)
        :param output_vcf_path: The path for the output VCF file (or "-" to write to stdout)
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
        :raise VcfAnnotatorAbortError: if a child class aborts annotation while a
            record is being processed
        """
        self.raise_for_output_args(output_vcf_path, **kwargs)
        # This can be pushed up to the click arg parsing too
        pysam_in_filename = (
            "-" if input_vcf_path == "-" else str(input_vcf_path.absolute())
        )
        vcf = pysam.VariantFile(filename=pysam_in_filename, mode="r")
        if output_vcf_path:
            self._update_vcf_header(vcf, compute_for_ref, vrs_attributes)
            pysam_out_filename = (
                "-" if output_vcf_path == "-" else str(output_vcf_path.absolute())
            )
            vcf_out = pysam.VariantFile(pysam_out_filename, mode="w", header=vcf.header)
        else:
            vcf_out = None

        allele_collection = (
            []
            if (self.collect_alleles and self.should_collect_alleles(**kwargs))
            else None
        )
        try:
            for record in vcf:
                if vcf_out:
                    additional_info_fields = [FieldName.IDS_FIELD]
                    if vrs_attributes:
                        additional_info_fields += [
                            FieldName.STARTS_FIELD,
                            FieldName.ENDS_FIELD,
                            FieldName.STATES_FIELD,
                            FieldName.LENGTHS_FIELD,
                            FieldName.REPEAT_SUBUNIT_LENGTHS_FIELD,
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
                except VcfAnnotatorAbortError:
                    _logger.exception(
                        "Aborting VCF annotation at %s-%s", record.chrom, record.pos
                    )
                    raise
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
                        # pysam correctly converts None → "." for Integer/Float
                        # INFO fields, but for String fields it converts None →
                        # "" (empty bytes), violating the VCF spec.  Work around
                        # by using the literal string "." for String-typed fields.
                        if k in _STRING_FIELDS:
                            record.info[k.value] = [
                                "." if v in ("", None) else v
                                for v in vrs_field_data[k.value]
                            ]
                        else:
                            record.info[k.value] = [
                                None if v in ("", None) else v
                                for v in vrs_field_data[k.value]
                            ]
                    vcf_out.write(record)
        finally:
            vcf.close()
            if vcf_out:
                vcf_out.close()

        if allele_collection is not None:
            self.on_vrs_object_collection(allele_collection, **kwargs)

    @abc.abstractmethod
    def on_vrs_object(
        self, vcf_coords: str, vrs_allele: Allele, **kwargs
    ) -> Allele | None:
        """Perform side-effects (eg additional annotation or storage) or additional
        filtering on VRS alleles as they are constructed during VCF annotation.

        Reimplement in a child class to add custom logic. Otherwise, simply pass through
        ``vrs_allele`` without altering it further or storing it.

        Exceptions raised here are recorded in the VCF record's ``VRS_Error`` INFO
        field and annotation continues with the next record, unless the exception is
        a ``VcfAnnotatorAbortError``, which stops annotation altogether.

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

        :param vrs_alleles_collection: VRS alleles constructed from ingested VCF
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
        :param allele_collection: All constructed VRS objects. Can be `None` if no data dumps
            will be created.
        :param vrs_field_data: If `vrs_data`, keys are VRS Fields and values are list
            of VRS data. Empty dict otherwise.
        :param assembly: The assembly used in `vcf_coords`
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
        if vrs_obj is not None:
            vrs_obj = self.on_vrs_object(vcf_coords, vrs_obj, **kwargs)
        else:
            _logger.debug(
                "None was returned when translating %s from gnomad", vcf_coords
            )

        if allele_collection is not None and vrs_obj:
            allele_collection.append(vrs_obj)

        if vrs_field_data:
            allele_id = vrs_obj.id if vrs_obj else ""
            vrs_field_data[FieldName.IDS_FIELD].append(allele_id)

            if vrs_attributes:
                # Initialize fields with None for missing values
                # pysam will convert None to "." in VCF output
                start = end = None
                alt = None
                length = repeat_subunit_length = None

                if vrs_obj:
                    # Common fields for all state types
                    start = vrs_obj.location.start
                    end = vrs_obj.location.end
                    state = vrs_obj.state
                    state_type = state.type

                    if state_type == VrsType.LIT_SEQ_EXPR:
                        # Sequence is required
                        alt = state.sequence.root
                    elif state_type == VrsType.REF_LEN_EXPR:
                        # Length is required, sequence is optional
                        length = state.length
                        if length is None:
                            err_msg = f"{state_type} requires a non-empty length: {vcf_coords}"
                            raise VcfAnnotatorError(err_msg)
                        if isinstance(length, Range):
                            err_msg = f"{state_type} with Range length not supported for VCF annotation: {vcf_coords}"
                            raise VcfAnnotatorError(err_msg)
                        repeat_subunit_length = state.repeatSubunitLength
                        # Only include sequence if within rle_seq_limit
                        if state.sequence is not None and (
                            RLE_SEQ_LIMIT is None or length <= RLE_SEQ_LIMIT
                        ):
                            alt = state.sequence.root
                    elif state_type == VrsType.LEN_EXPR:
                        # Length is required, no sequence field
                        length = state.length
                        if length is None:
                            err_msg = f"{state_type} requires a non-empty length: {vcf_coords}"
                            raise VcfAnnotatorError(err_msg)
                        if isinstance(length, Range):
                            err_msg = f"{state_type} with Range length not supported for VCF annotation: {vcf_coords}"
                            raise VcfAnnotatorError(err_msg)
                    else:
                        err_msg = f"Unsupported state type '{state_type}' for VCF annotation: {vcf_coords}"
                        raise VcfAnnotatorError(err_msg)

                vrs_field_data[FieldName.STARTS_FIELD].append(start)
                vrs_field_data[FieldName.ENDS_FIELD].append(end)
                vrs_field_data[FieldName.STATES_FIELD].append(alt)
                vrs_field_data[FieldName.LENGTHS_FIELD].append(length)
                vrs_field_data[FieldName.REPEAT_SUBUNIT_LENGTHS_FIELD].append(
                    repeat_subunit_length
                )

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
        :param allele_collection: List containing the VRS object information for
            the VCF, if `self.collect_alleles` is `True`.
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

    collect_alleles = True
    pkl_arg_name = "output_pkl_path"
    ndjson_arg_name = "output_ndjson_path"

    def should_collect_alleles(self, **kwargs) -> bool:
        """Inhibit allele collection parameter if no means of output are given.

        :kwparam output_pkl_path: Optional path to output PKL dump of all alleles
        :kwparam output_ndjson_path: Optional path to output NDJSON dump of all alleles
        :return: ``True`` if at least one of the output path args is ``True``
        """
        return bool(kwargs.get("output_pkl_path") or kwargs.get("output_ndjson_path"))

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
        a PKL or NDJSON dump of the corresponding VRS objects.

        :param input_vcf_path: Location of input VCF
        :param output_vcf_path: The path for the output VCF file
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
        :kwparam output_pkl_path: Optional path to output PKL dump of all alleles
        :kwparam output_ndjson_path: Optional path to output NDJSON dump of all alleles
        :raise VCFAnnotatorError: if no output formats are selected
        """
        return super().annotate(
            input_vcf_path,
            output_vcf_path,
            vrs_attributes,
            assembly,
            compute_for_ref,
            require_validation,
            **kwargs,
        )

    def raise_for_output_args(self, output_vcf_path: Path | None, **kwargs) -> None:
        """Raise an exception if no output (VCF, PKL, NDJSON) appears to be configured
        or declared.

        :param output_vcf_path: VCF output path arg passed to `annotate()`
        :kwparam output_pkl_path: optional path to PKL output
        :kwparam output_ndjson_path: optional path to NDJSON output
        :raise VCFAnnotatorArgsError: if no output args are given
        """
        if (
            output_vcf_path is None
            and kwargs.get(self.pkl_arg_name) is None
            and kwargs.get(self.ndjson_arg_name) is None
        ):
            msg = f"No VCF, PKL, or NDJSON output path provided -- must pass at least one of `output_vcf_path`, `{self.pkl_arg_name}`, `{self.ndjson_arg_name}` to annotate()."
            raise VcfAnnotatorArgsError(msg)

    def on_vrs_object(
        self,
        vcf_coords: str,  # noqa: ARG002
        vrs_allele: Allele,
        **kwargs,  # noqa: ARG002
    ) -> Allele | None:
        """Perform side-effects (eg additional annotation or storage) or additional
        filtering on VRS alleles as they are constructed during VCF annotation.

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

        :param vrs_alleles_collection: VRS alleles constructed from ingested VCF
        """
        if vrs_alleles_collection is not None:
            output_pkl_path = kwargs.get(self.pkl_arg_name)
            if output_pkl_path is not None:
                dump_alleles_to_pkl(vrs_alleles_collection, output_pkl_path)

            output_ndjson_path = kwargs.get(self.ndjson_arg_name)
            if output_ndjson_path is not None:
                dump_alleles_to_ndjson(vrs_alleles_collection, output_ndjson_path)


class VcfReader:
    """Read VRS Alleles back from a VCF annotated with VRS attributes.

    The input VCF must have been annotated with VRS attributes (e.g. by
    ``VcfAnnotator.annotate()`` with ``vrs_attributes=True``, or by
    ``vrs-annotate vcf --vrs-attributes``). Each VRS Allele is reconstructed from the
    record's CHROM and its VRS INFO fields, and its computed identifier must match the
    annotated ``VRS_Allele_IDs`` value.
    """

    def __init__(self, data_proxy: _DataProxy) -> None:
        """Initialize the VcfReader class.

        :param data_proxy: GA4GH sequence dataproxy instance, used to look up the
            refget accession for each CHROM
        """
        self.data_proxy = data_proxy

    def read(
        self, input_vcf_path: Path | Literal["-"], assembly: str = "GRCh38"
    ) -> Iterator[tuple[pysam.VariantRecord, list[Allele | None]]]:
        """Iterate over the records of an annotated VCF along with their VRS Alleles.

        For each record, yield the pysam record and a list of identified VRS Alleles,
        one per ``VRS_Allele_IDs`` value: one per REF and ALT allele
        (``record.alleles``) if the VCF was annotated with REF alleles, or one per ALT
        allele (``record.alts``) if not. The list holds
        ``None`` for alleles without a VRS ID, such as ``*`` ALTs, and for every
        allele of a record that failed annotation (see its ``VRS_Error`` INFO field).

        The input VCF is opened when iteration starts and closed when iteration ends.

        ReferenceLengthExpression sequences longer than the annotator's sequence
        limit aren't stored in the VCF, so these states are returned without a
        ``sequence``. This doesn't change their VRS identifiers.

        :param input_vcf_path: Location of input VCF (or "-" to read from stdin)
        :param assembly: The assembly used in `input_vcf_path` data
        :raise VcfReaderError: if the VCF wasn't annotated with VRS attributes, or if
            the VRS Alleles of a record can't be reconstructed from its annotations,
            including when a computed identifier doesn't match the annotated one
        """
        pysam_in_filename = (
            "-" if input_vcf_path == "-" else str(input_vcf_path.absolute())
        )
        with pysam.VariantFile(filename=pysam_in_filename, mode="r") as vcf:
            missing_fields = [
                field.value
                for field in (FieldName.IDS_FIELD, *_VRS_ATTRIBUTE_FIELDS)
                if field.value not in vcf.header.info
            ]
            if missing_fields:
                msg = f"VCF is missing VRS annotation INFO fields {', '.join(missing_fields)} -- it must be annotated with VRS attributes (`--vrs-attributes`) to read VRS Alleles from it."
                raise VcfReaderError(msg)
            includes_ref = vcf.header.info[FieldName.IDS_FIELD.value].number == "R"

            for record in vcf:
                yield record, self._get_vrs_alleles(record, assembly, includes_ref)

    def _get_vrs_alleles(
        self, record: pysam.VariantRecord, assembly: str, includes_ref: bool
    ) -> list[Allele | None]:
        """Reconstruct the VRS Alleles for a record from its VRS annotations.

        :param record: A row in the VCF file
        :param assembly: The assembly used in `record`
        :param includes_ref: Whether the VCF was annotated with VRS Alleles for REFs
        :return: VRS Alleles, or ``None`` for alleles without a VRS ID, in the order of
            ``record.alleles`` if `includes_ref`, otherwise of ``record.alts``
        :raise VcfReaderError: if the VRS annotations of `record` are incomplete or
            inconsistent
        """
        num_alleles = len(record.alleles) if includes_ref else len(record.alts or ())
        if FieldName.IDS_FIELD.value not in record.info:
            return [None] * num_alleles

        vcf_coords = f"{record.chrom}-{record.pos}"
        ids, starts, ends, states, lengths, repeat_subunit_lengths = (
            self._get_info_values(record, field, num_alleles, vcf_coords)
            for field in (FieldName.IDS_FIELD, *_VRS_ATTRIBUTE_FIELDS)
        )
        refget_accession = None
        vrs_alleles = []
        for vrs_id, start, end, sequence, length, repeat_subunit_length in zip(
            ids, starts, ends, states, lengths, repeat_subunit_lengths, strict=True
        ):
            if vrs_id is None:
                vrs_alleles.append(None)
                continue

            if start is None or end is None:
                msg = f"VRS Allele {vrs_id} at {vcf_coords} is missing its {FieldName.STARTS_FIELD.value} or {FieldName.ENDS_FIELD.value} value"
                raise VcfReaderError(msg)

            if refget_accession is None:
                sequence_id = f"{assembly}:{record.chrom}"
                refget_accession = self.data_proxy.derive_refget_accession(sequence_id)
                if refget_accession is None:
                    msg = f"Unable to find refget accession for {sequence_id}"
                    raise VcfReaderError(msg)

            try:
                vrs_allele = Allele(
                    location=SequenceLocation(
                        sequenceReference=SequenceReference(
                            refgetAccession=refget_accession
                        ),
                        start=start,
                        end=end,
                    ),
                    state=self._get_state(sequence, length, repeat_subunit_length),
                )
            except ValidationError as e:
                msg = f"Invalid VRS annotations for VRS Allele {vrs_id} at {vcf_coords}: {e}"
                raise VcfReaderError(msg) from e
            vrs_allele.id = ga4gh_identify(vrs_allele)
            vrs_allele.location.id = ga4gh_identify(vrs_allele.location)
            if vrs_allele.id != vrs_id:
                msg = f"VRS Allele reconstructed at {vcf_coords} has ID {vrs_allele.id}, but was annotated with {vrs_id}"
                raise VcfReaderError(msg)
            vrs_alleles.append(vrs_allele)
        return vrs_alleles

    @staticmethod
    def _get_info_values(
        record: pysam.VariantRecord, field: FieldName, num_alleles: int, vcf_coords: str
    ) -> tuple:
        """Get the per-allele values of a VRS INFO field, with missing values as
        ``None``.

        :param record: A row in the VCF file
        :param field: VRS INFO field to get
        :param num_alleles: Number of values the field must have
        :param vcf_coords: CHR-POS of `record`, for error messages
        :return: The values of `field`, or ``None`` for each allele if `record` doesn't
            have it
        :raise VcfReaderError: if `field` doesn't have `num_alleles` values
        """
        values = record.info.get(field.value)
        if values is None:
            return (None,) * num_alleles
        if len(values) != num_alleles:
            msg = f"Expected {num_alleles} {field.value} values at {vcf_coords}, but found {len(values)}"
            raise VcfReaderError(msg)
        # Missing String values are written as "." (or "" by older VRS-Python versions)
        return tuple(None if value in ("", ".") else value for value in values)

    @staticmethod
    def _get_state(
        sequence: str | None, length: int | None, repeat_subunit_length: int | None
    ) -> LiteralSequenceExpression | ReferenceLengthExpression | LengthExpression:
        """Construct a VRS Allele state from its VRS annotation values.

        :param sequence: ``VRS_States`` value
        :param length: ``VRS_Lengths`` value
        :param repeat_subunit_length: ``VRS_RepeatSubunitLengths`` value
        :return: LiteralSequenceExpression if there's no `length`, otherwise
            ReferenceLengthExpression if there's a `repeat_subunit_length`, otherwise
            LengthExpression
        """
        # An empty sequence is written as a missing value, which is also used for
        # ReferenceLengthExpression sequences left out because of their length.
        if length is None:
            return LiteralSequenceExpression(sequence=sequence or "")
        if repeat_subunit_length is None:
            return LengthExpression(length=length)
        if sequence is None and length == 0:
            sequence = ""
        return ReferenceLengthExpression(
            length=length, repeatSubunitLength=repeat_subunit_length, sequence=sequence
        )
