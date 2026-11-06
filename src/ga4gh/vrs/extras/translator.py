"""Translates various external formats into VRS models.

Input formats: VRS (serialized), hgvs, spdi, gnomad (vcf)
Output formats: VRS (serialized), hgvs, spdi, gnomad (vcf)

"""

import logging
import re
from abc import ABC
from collections.abc import Mapping
from enum import Enum
from typing import Protocol

from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models, normalize
from ga4gh.vrs.dataproxy import SequenceProxy, _DataProxy
from ga4gh.vrs.extras.decorators import lazy_property
from ga4gh.vrs.normalize import denormalize_reference_length_expression
from ga4gh.vrs.utils.hgvs_tools import HgvsTools

_logger = logging.getLogger(__name__)


class TranslationDefault(Enum):
    """Mark options that use the translator instance setting.

    USE_INSTANCE_DEFAULT distinguishes an omitted option from an explicit None.
    For rle_seq_limit, the default uses the instance limit, while None means no
    limit. Enum members preserve their identity through copying and pickling.
    """

    INSTANCE = "instance"


USE_INSTANCE_DEFAULT = TranslationDefault.INSTANCE
"""Use the translator instance setting for an omitted option."""


class VariationToStrProtocol(Protocol):
    """Protocol for translating VRS objects to other string expressions.

    This protocol defines a callable interface for translating a VRS object
    into variation strings, with optional keyword arguments for customization.
    """

    def __call__(self, vo: models._VariationBase, **kwargs) -> list[str]:
        """Translate vrs object `vo` to variation string expressions"""


class VariationFromStrProtocol(Protocol):
    """Protocol for translating variation strings to VRS objects.

    This protocol defines a callable interface for translating a variation
    string into a VRS object, with optional keyword arguments for customization.
    """

    def __call__(self, expr: str, **kwargs) -> models._VariationBase | None:
        """Translate variation string `expr` to a VRS object"""


class _Translator(ABC):  # noqa: B024
    """abstract class / interface for VRS to/from translation needs

     Translates various variation formats to and from GA4GH VRS models

    Format-specific `from_*` parsers follow this pattern:
    * If the argument does not appear to be an appropriate type, None is returned
    * Otherwise, the argument is expected to be of the correct type.  If an error occurs during processing,
      an exception is raised.
    * Otherwise, the VRS object is returned

    """

    gnomad_re = re.compile(
        r"(?P<chr>[^-]+)-(?P<pos>\d+)-(?P<ref>[ACGTURYKMSWBDHVN]+)-(?P<alt>[ACGTURYKMSWBDHVN]+)",
        re.IGNORECASE,
    )
    spdi_re = re.compile(
        r"(?P<ac>[^:]+):(?P<pos>\d+):(?P<del_len_or_seq>\w*):(?P<ins_seq>\w*)"
    )

    def __init__(
        self,
        data_proxy: _DataProxy,
        default_assembly_name: str = "GRCh38",
        identify: bool = True,
        rle_seq_limit: int | None = 50,
    ) -> None:
        """Initialize shared translation settings and format dispatch dictionaries.

        :param data_proxy: Sequence data provider used during translation.
        :param default_assembly_name: Default assembly for inputs without an assembly.
        :param identify: Whether to assign GA4GH identifiers to translated objects.
        :param rle_seq_limit: Maximum RLE sequence length to include, or None for no limit.
        """
        self.default_assembly_name = default_assembly_name
        self.data_proxy = data_proxy
        self.identify = identify
        self.rle_seq_limit = rle_seq_limit
        self.to_vrs_translators: dict[str, VariationFromStrProtocol] = {}
        self.from_vrs_translators: dict[str, VariationToStrProtocol] = {}

    def to_vrs(
        self, var: str | dict, fmt: str | None = None, **kwargs
    ) -> models._VariationBase:
        """Translate a variation expression or VRS JSON dictionary to a VRS object.

        If fmt is omitted, try the registered parsers in order. Format-specific
        methods return None for unrecognized input; this method raises ValueError.

        :param var: Variation expression, or a VRS JSON dictionary with fmt="vrs".
        :param fmt: Input format, or None to try the registered parsers.
        :param kwargs: Options forwarded to the selected parser. Unused options
            are ignored.

            For AlleleTranslator HGVS, SPDI, and gnomAD inputs:

            * ``do_normalize`` (bool): Whether to perform fully justified
              normalization. Defaults to True.
            * ``rle_seq_limit`` (int or None): Maximum sequence length to include
              when normalization produces a ReferenceLengthExpression. Omission
              or ``USE_INSTANCE_DEFAULT`` uses the instance's ``rle_seq_limit``.
              Set to 0 to exclude the sequence, or None for no limit.

            Additional options for AlleleTranslator gnomAD inputs:

            * ``assembly_name`` (str or None): Assembly used for the input.
              Omission or None uses the instance's ``default_assembly_name``.
            * ``require_validation`` (bool): Whether reference mismatches raise
              DataProxyValidationError. Defaults to True. False allows reference
              mismatches; location bounds are still validated.

            For CnvTranslator HGVS inputs:

            * ``copies`` (int, models.Range, or None): Copy count or range,
              including zero. If supplied, returns CopyNumberCount and takes
              precedence over ``copy_change``. Defaults to None.
            * ``copy_change`` (models.CopyChange or None): Copy change used when
              ``copies`` is omitted or None. Defaults to LOSS (EFO:0030067) for
              deletions and GAIN (EFO:0030070) for duplications.
        :return: Translated VRS variation.
        :raises NotImplementedError: If the requested format is unsupported.
        :raises ValueError: If the input cannot be parsed or represented.
        :raises DataProxyValidationError: If required reference or location validation fails.
        """
        if fmt:
            try:
                t = self.to_vrs_translators[fmt]
            except KeyError as e:
                msg = f"{fmt} is not supported"
                raise NotImplementedError(msg) from e
            else:
                o = t(var, **kwargs)
                if o is None:
                    msg = f"Unable to parse data as {fmt} variation"
                    raise ValueError(msg)
                return o

        for t in self.to_vrs_translators.values():
            o = t(var, **kwargs)
            if o:
                return o

        formats = list(self.to_vrs_translators.keys())
        msg = f"Unable to parse data as {', '.join(formats)}"
        raise ValueError(msg)

    def from_vrs(self, vo: models._VariationBase, fmt: str, **kwargs) -> list[str]:
        """Translate a VRS object to expressions in the requested format.

        :param vo: VRS variation to translate.
        :param fmt: Output format.
        :param kwargs: Options forwarded to the selected exporter. Unused options
            are ignored.

            For AlleleTranslator HGVS and SPDI exports:

            * ``namespace`` (str or None): Sequence identifier namespace.
              Defaults to "refseq". None includes all available aliases.

            Additional options for AlleleTranslator SPDI exports:

            * ``ref_seq_limit`` (int or None): Maximum reference sequence length
              to include in the SPDI expression. Longer reference sequences are
              represented by their length. Defaults to 0; None always includes
              the reference sequence. The alternate sequence is always included,
              including when derived from a ReferenceLengthExpression.
        :return: Expressions in the requested format.
        :raises NotImplementedError: If the requested format is unsupported.
        :raises ValueError: If the object cannot be represented in the requested format.
        """
        try:
            t = self.from_vrs_translators[fmt]
        except KeyError as e:
            msg = f"{fmt} is not supported"
            raise NotImplementedError(msg) from e

        return t(vo, **kwargs)

    def from_vrs_dict(self, var: dict, **kwargs) -> models._VariationBase | None:  # noqa: ARG002
        """Parse a VRS JSON dictionary into a VRS variation object.

        :param var: Dictionary containing a VRS variation and its ``type``.
        :param kwargs: Unused options accepted for format dispatch.
        :return: Parsed variation, or None for non-mappings or missing ``type``.
        :raises pydantic.ValidationError: If the dictionary is not a valid variation.
        """
        if not isinstance(var, Mapping):
            return None
        if "type" not in var:
            return None
        return models.Variation.model_validate(var).root

    ############################################################################
    # INTERNAL

    @lazy_property
    def hgvs_tools(self) -> HgvsTools:
        """Instantiate and return an HgvsTools instance"""
        return HgvsTools(self.data_proxy)


class AlleleTranslator(_Translator):
    """Class for translating formats to and from VRS Alleles"""

    def __init__(
        self,
        data_proxy: _DataProxy,
        default_assembly_name: str = "GRCh38",
        identify: bool = True,
    ) -> None:
        """Initialize translation to and from VRS Alleles.

        :param data_proxy: Sequence data provider used during translation.
        :param default_assembly_name: Default assembly for gnomAD inputs.
        :param identify: Whether to assign GA4GH identifiers to translated objects.
        """
        super().__init__(data_proxy, default_assembly_name, identify)

        self.to_vrs_translators = {
            "gnomad": self.from_gnomad,
            "hgvs": self.from_hgvs,
            "spdi": self.from_spdi,
            "vrs": self.from_vrs_dict,
        }

        self.from_vrs_translators = {
            "hgvs": self.to_hgvs,
            "spdi": self.to_spdi,
        }

    def _create_allele(
        self,
        values: dict,
        *,
        do_normalize: bool = True,
        rle_seq_limit: int | None | TranslationDefault = USE_INSTANCE_DEFAULT,
    ) -> models.Allele:
        """Build and post-process an Allele from parsed coordinates and sequence.

        :param values: Parsed ``sequence_id``, ``refget_accession``, ``start``,
            ``end``, and ``literal_sequence`` values.
        :param do_normalize: Whether to perform fully justified normalization (default True).
        :param rle_seq_limit: Maximum RLE sequence length to include.
            USE_INSTANCE_DEFAULT uses the instance limit; 0 excludes the sequence
            and None includes it without a limit.
        :return: Constructed Allele, normalized and identified as configured.
        :raises DataProxyValidationError: If the location is outside the reference sequence.
        """
        self.data_proxy.validate_location_bounds(
            values["sequence_id"], values["start"], values["end"]
        )
        seq_ref = models.SequenceReference(refgetAccession=values["refget_accession"])
        location = models.SequenceLocation(
            sequenceReference=seq_ref, start=values["start"], end=values["end"]
        )
        state = models.LiteralSequenceExpression(sequence=values["literal_sequence"])
        allele = models.Allele(location=location, state=state)
        return self._post_process_imported_allele(
            allele, do_normalize=do_normalize, rle_seq_limit=rle_seq_limit
        )

    def from_gnomad(
        self,
        gnomad_expr: str,
        *,
        assembly_name: str | None = None,
        require_validation: bool = True,
        do_normalize: bool = True,
        rle_seq_limit: int | None | TranslationDefault = USE_INSTANCE_DEFAULT,
        **kwargs,  # noqa: ARG002
    ) -> models.Allele | None:
        """Parse gnomAD-style VCF expression into VRS Allele

        :param gnomad_expr: Variant in chr-pos-ref-alt format.
        :param assembly_name: Assembly used for the input, or None for the instance default.
        :param require_validation: Whether reference mismatches raise an error (default True).
        :param do_normalize: Whether to perform fully justified normalization (default True).
        :param rle_seq_limit: Maximum RLE sequence length to include.
            USE_INSTANCE_DEFAULT uses the instance limit; 0 excludes the sequence
            and None includes it without a limit.
        :param kwargs: Unused options accepted for format dispatch.
        :return: Translated Allele, or None if the expression or accession is unrecognized.
        :raises DataProxyValidationError: If reference validation is required and fails,
            or the location is outside the reference sequence.

        #>>> a = tlr.from_gnomad("1-55516888-G-GA")
        #>>> a.model_dump()
        {
          'location': {
            'end': 55516888,
            'start': 55516887,
            'sequenceReference': {
              'refgetAccession': 'SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO',
              'type': 'SequenceReference'
            },
            'type': 'SequenceLocation'
          },
          'state': {
            'sequence': 'GA',
            'type': 'LiteralSequenceExpression'
          },
          'type': 'Allele'
        }

        """
        if not isinstance(gnomad_expr, str):
            return None
        m = self.gnomad_re.fullmatch(gnomad_expr)
        if not m:
            return None

        g = m.groupdict()
        if assembly_name is None:
            assembly_name = self.default_assembly_name
        sequence = assembly_name + ":" + g["chr"]
        refget_accession = self.data_proxy.derive_refget_accession(sequence)
        if not refget_accession:
            return None

        start = int(g["pos"]) - 1
        ref = g["ref"].upper()
        alt = g["alt"].upper()
        end = start + len(ref)
        ins_seq = alt

        # validation checks
        self.data_proxy.validate_ref_seq(
            sequence,
            start,
            end,
            ref,
            require_validation=require_validation,
        )

        values = {
            "sequence_id": sequence,
            "refget_accession": refget_accession,
            "start": start,
            "end": end,
            "literal_sequence": ins_seq,
        }
        return self._create_allele(
            values, do_normalize=do_normalize, rle_seq_limit=rle_seq_limit
        )

    def from_hgvs(
        self,
        hgvs_expr: str,
        *,
        do_normalize: bool = True,
        rle_seq_limit: int | None | TranslationDefault = USE_INSTANCE_DEFAULT,
        **kwargs,  # noqa: ARG002
    ) -> models.Allele | None:
        """Parse an HGVS expression into a VRS Allele.

        :param hgvs_expr: HGVS variant expression.
        :param do_normalize: Whether to perform fully justified normalization (default True).
        :param rle_seq_limit: Maximum RLE sequence length to include.
            USE_INSTANCE_DEFAULT uses the instance limit; 0 excludes the sequence
            and None includes it without a limit.
        :param kwargs: Unused options accepted for format dispatch.
        :return: Translated Allele, or None if the expression or accession is unrecognized.
        :raises ValueError: If HGVS positions or edits cannot be represented as an Allele.
        :raises DataProxyValidationError: If the location is outside the reference sequence.
        """
        sv = self.hgvs_tools.parse(hgvs_expr)
        if sv is not None and self.hgvs_tools.has_uncertain_range(sv):
            msg = (
                "Uncertain-range HGVS expressions are not supported for Allele "
                "translation; use CnvTranslator for del/dup"
            )
            raise ValueError(msg)

        allele_values = self.hgvs_tools.extract_allele_values(hgvs_expr)
        if allele_values:
            return self._create_allele(
                allele_values, do_normalize=do_normalize, rle_seq_limit=rle_seq_limit
            )
        return None

    def from_spdi(
        self,
        spdi_expr: str,
        *,
        do_normalize: bool = True,
        rle_seq_limit: int | None | TranslationDefault = USE_INSTANCE_DEFAULT,
        **kwargs,  # noqa: ARG002
    ) -> models.Allele | None:
        """Parse SPDI expression in to a GA4GH Allele

        :param spdi_expr: Variant in accession:position:deletion:insertion format.
        :param do_normalize: Whether to perform fully justified normalization (default True).
        :param rle_seq_limit: Maximum RLE sequence length to include.
            USE_INSTANCE_DEFAULT uses the instance limit; 0 excludes the sequence
            and None includes it without a limit.
        :param kwargs: Unused options accepted for format dispatch.
        :return: Translated Allele, or None if the expression or accession is unrecognized.
        :raises DataProxyValidationError: If the location is outside the reference sequence.

        #>>> a = tlr.from_spdi("NC_000013.11:32936731:1:C")
        #>>> a.model_dump()
        {
          'location': {
            'end': 32936732,
            'start': 32936731,
            'sequenceReference': {
                'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT',
                'type': 'SequenceReference'
            },
            'type': 'SequenceLocation'
          },
          'state': {
              'sequence': 'C',
              'type': 'LiteralSequenceExpression'
          },
          'type': 'Allele'
        }
        """
        if not isinstance(spdi_expr, str):
            return None

        m = self.spdi_re.fullmatch(spdi_expr)
        if not m:
            return None

        g = m.groupdict()
        refget_accession = self.data_proxy.derive_refget_accession(g["ac"])
        if not refget_accession:
            return None

        start = int(g["pos"])
        try:
            del_len = int(g["del_len_or_seq"])
        except ValueError:
            del_len = len(g["del_len_or_seq"])
        end = start + del_len
        ins_seq = g["ins_seq"]

        values = {
            "sequence_id": g["ac"],
            "refget_accession": refget_accession,
            "start": start,
            "end": end,
            "literal_sequence": ins_seq,
        }

        return self._create_allele(
            values, do_normalize=do_normalize, rle_seq_limit=rle_seq_limit
        )

    def to_hgvs(
        self,
        vo: models.Allele,
        namespace: str | None = "refseq",
        **kwargs,  # noqa: ARG002
    ) -> list[str]:
        """Generate HGVS expressions for a VRS Allele.

        :param vo: Allele to translate.
        :param namespace: Sequence identifier namespace, or None for all aliases.
        :param kwargs: Unused options accepted for format dispatch.
        :return: HGVS expressions, or an empty list if no aliases are available.
        :raises ValueError: If the object cannot be represented as HGVS.
        """
        return self.hgvs_tools.from_allele(vo, namespace)

    def to_spdi(
        self,
        vo: models.Allele,
        namespace: str | None = "refseq",
        *,
        ref_seq_limit: int | None = 0,
        **kwargs,  # noqa: ARG002
    ) -> list[str]:
        """Generate a *list* of SPDI expressions for VRS Allele.

        If `namespace` is not None, returns SPDI strings for the
        specified namespace.

        If `namespace` is None, returns SPDI strings for all alias
        translations.

        If `ref_seq_limit` is specified, the reference sequence is
        included in the SPDI expression only if it is below the limit.
        Otherwise only the length of the reference sequence is
        included. If the limit is None, the reference sequence is
        always included. In all cases, the alt sequence is included.
        Default is 0 (never include reference sequence).

        If no alias translations are available, an empty list is
        returned.

        If the VRS object cannot be expressed as SPDI, raises ValueError.

        SPDI and VRS use identical normalization. The incoming Allele
        is expected to be normalized per VRS spec.

        :param vo: Allele to translate.
        :param namespace: Sequence identifier namespace, or None for all aliases.
        :param ref_seq_limit: Maximum reference sequence length to include (default 0).
            None always includes the reference sequence.
        :param kwargs: Unused options accepted for format dispatch.
        :return: SPDI expressions, or an empty list if no aliases are available.
        :raises ValueError: If the object cannot be represented as SPDI.
        """
        sequence = f"ga4gh:{vo.location.get_refget_accession()}"
        aliases = self.data_proxy.translate_sequence_identifier(sequence, namespace)
        aliases = [a.split(":")[1] for a in aliases]
        seq_proxies = {a: SequenceProxy(self.data_proxy, a) for a in aliases}
        start, end = vo.location.start, vo.location.end
        spdi_exprs = []

        for alias in aliases:
            # Get the reference sequence
            seq_proxy = seq_proxies[alias]
            ref_seq = seq_proxy[start:end]

            if vo.state.type == models.VrsType.REF_LEN_EXPR.value:
                # Derived from reference. sequence included if under limit, but
                # we can derive it again from the reference.
                alt_seq = denormalize_reference_length_expression(
                    ref_seq=ref_seq,
                    repeat_subunit_length=vo.state.repeatSubunitLength,
                    alt_length=vo.state.length,
                )
                # Warn if the derived sequence is different from the one in the object
                if vo.state.sequence and vo.state.sequence.root != alt_seq:
                    _logger.warning(
                        "Derived sequence '%s' is different from provided state.sequence '%s'",
                        alt_seq,
                        vo.state.sequence.root,
                    )
            else:
                alt_seq = vo.state.sequence.root

            # Optionally allow using the length of the reference sequence
            # instead of the sequence itself.
            if ref_seq_limit is not None and len(ref_seq) > int(ref_seq_limit):
                ref_seq = len(ref_seq)

            spdi_expr = f"{alias}:{start}:{ref_seq}:{alt_seq}"
            spdi_exprs.append(spdi_expr)

        return spdi_exprs

    def _post_process_imported_allele(
        self,
        allele: models.Allele,
        *,
        do_normalize: bool = True,
        rle_seq_limit: int | None | TranslationDefault = USE_INSTANCE_DEFAULT,
    ) -> models.Allele:
        """Normalize an imported Allele and assign identifiers as configured.

        :param allele: Allele to process. Identifiers are assigned in place.
        :param do_normalize: Whether to perform fully justified normalization (default True).
        :param rle_seq_limit: Maximum RLE sequence length to include.
            USE_INSTANCE_DEFAULT uses the instance limit; 0 excludes the sequence
            and None includes it without a limit.
        :return: Processed Allele. Normalization may return a new object.
        """
        if rle_seq_limit is USE_INSTANCE_DEFAULT:
            rle_seq_limit = self.rle_seq_limit

        if do_normalize:
            allele = normalize(
                allele,
                self.data_proxy,
                rle_seq_limit=rle_seq_limit,
            )

        if self.identify:
            allele.id = ga4gh_identify(allele)
            allele.location.id = ga4gh_identify(allele.location)

        return allele


class CnvTranslator(_Translator):
    """Class for translating formats from format to VRS Copy Number"""

    def __init__(
        self,
        data_proxy: _DataProxy,
        default_assembly_name: str = "GRCh38",
        identify: bool = True,
    ) -> None:
        """Initialize translation from HGVS to VRS copy number objects.

        :param data_proxy: Sequence data provider used during translation.
        :param default_assembly_name: Default assembly name.
        :param identify: Whether to assign GA4GH identifiers to translated objects.
        """
        super().__init__(data_proxy, default_assembly_name, identify)
        self.to_vrs_translators = {
            "hgvs": self.from_hgvs,
        }

    def from_hgvs(
        self,
        hgvs_dup_del_expr: str,
        *,
        copies: int | models.Range | None = None,
        copy_change: models.CopyChange | None = None,
        **kwargs,  # noqa: ARG002
    ) -> models.CopyNumberChange | models.CopyNumberCount | None:
        """Parse an HGVS deletion or duplication into a VRS copy number object.

        :param hgvs_dup_del_expr: HGVS deletion or duplication expression.
        :param copies: Copy count or range, including zero. If supplied, returns
            CopyNumberCount and takes precedence over ``copy_change``.
        :param copy_change: Copy change, defaulting to LOSS for deletions and
            GAIN for duplications when ``copies`` is omitted.
        :param kwargs: Unused options accepted for format dispatch.
        :return: CopyNumberCount if copies is supplied, otherwise CopyNumberChange,
            or None if the expression or accession is unrecognized.
        :raises ValueError: If the edit is not a deletion or duplication, or its
            positions cannot be represented as a copy number location.
        :raises DataProxyValidationError: If the location is outside the reference sequence.
        """
        sv = self.hgvs_tools.parse(hgvs_dup_del_expr)
        if not sv:
            return None

        sv_type = self.hgvs_tools.get_edit_type(sv)
        if sv_type not in {"del", "dup"}:
            msg = "Must provide a 'del' or 'dup'"
            raise ValueError(msg)

        if self.hgvs_tools.is_intronic(sv):
            msg = "Intronic HGVS variants are not supported"
            raise ValueError(msg)

        refget_accession = self.data_proxy.derive_refget_accession(sv.ac)
        if not refget_accession:
            return None

        location = self.hgvs_tools.build_cnv_location(sv, refget_accession)
        self.data_proxy.validate_location_bounds(sv.ac, location.start, location.end)

        if copies is not None:
            cnv = models.CopyNumberCount(location=location, copies=copies)
        else:
            if not copy_change:
                copy_change = (
                    models.CopyChange.LOSS
                    if sv_type == "del"
                    else models.CopyChange.GAIN
                )
            cnv = models.CopyNumberChange(
                location=location,
                copyChange=copy_change,
            )

        return self._post_process_imported_cnv(cnv)

    def _post_process_imported_cnv(
        self, copy_number: models.CopyNumberChange | models.CopyNumberCount
    ) -> models.CopyNumberChange | models.CopyNumberCount:
        """Provide common post-processing for imported Copy Numbers IN-PLACE."""
        if self.identify:
            copy_number.id = ga4gh_identify(copy_number)
            copy_number.location.id = ga4gh_identify(copy_number.location)

        return copy_number
