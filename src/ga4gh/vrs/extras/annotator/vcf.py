"""Annotate VCFs with VRS."""

import logging
import pickle
from pathlib import Path
from typing import ClassVar

import pysam

from ga4gh.core import VrsObjectIdentifierIs, use_ga4gh_compute_identifier_when
from ga4gh.vrs.dataproxy import _DataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator

_logger = logging.getLogger(__name__)


class VCFAnnotatorError(Exception):
    """Raise for errors specific to the VCF annotation process"""


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
    VCF_ESCAPE_MAP: ClassVar = str.maketrans(
        {
            "%": "%25",
            ";": "%3B",
            ",": "%2C",
            "\r": "%0D",
            "\n": "%0A",
            "\t": "%09",
        }
    )

    def __init__(self, data_proxy: _DataProxy) -> None:
        """Initialize the VCFAnnotator class.

        :param data_proxy: GA4GH sequence dataproxy instance.
        """
        self.data_proxy = data_proxy
        self.tlr = AlleleTranslator(self.data_proxy)

    def _update_vcf_header(
        self,
        vcf: pysam.VariantFile,
        incl_ref_allele: bool,
        incl_vrs_attrs: bool,
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
            self.VRS_ALLELE_IDS_FIELD,
            info_field_num,
            "String",
            f"The computed identifiers for the GA4GH VRS Alleles corresponding to the GT indexes of the {info_field_desc} alleles",
        )
        vcf.header.info.add(
            self.VRS_ERROR_FIELD,
            ".",
            "String",
            "If an error occurred computing a VRS Identifier, the error message",
        )

        if incl_vrs_attrs:
            vcf.header.info.add(
                self.VRS_STARTS_FIELD,
                info_field_num,
                "String",
                f"Interresidue coordinates used as the location starts for the GA4GH VRS Alleles corresponding to the GT indexes of the {info_field_desc} alleles",
            )
            vcf.header.info.add(
                self.VRS_ENDS_FIELD,
                info_field_num,
                "String",
                f"Interresidue coordinates used as the location ends for the GA4GH VRS Alleles corresponding to the GT indexes of the {info_field_desc} alleles",
            )
            vcf.header.info.add(
                self.VRS_STATES_FIELD,
                info_field_num,
                "String",
                f"The literal sequence states used for the GA4GH VRS Alleles corresponding to the GT indexes of the {info_field_desc} alleles",
            )

    def _process_allele(
        self,
        vcf_coords: str,
        vrs_data: dict,
        annotations: dict,
        assembly: str,
        vrs_data_key: str | None = None,
        create_pickle: bool = True,
        incl_vrs_attrs: bool = False,
        require_validation: bool = True,
    ) -> None:
        """Get VRS object given `vcf_coords`. `vrs_data` and `vrs_field_data` will
        be mutated.

        # TODO update
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

        if create_pickle and vrs_obj:
            key = vrs_data_key if vrs_data_key else vcf_coords
            vrs_data[key] = str(vrs_obj.model_dump(exclude_none=True))

        if annotations:
            allele_id = vrs_obj.id if vrs_obj else ""
            annotations[self.VRS_ALLELE_IDS_FIELD].append(allele_id)

            if incl_vrs_attrs:
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

                annotations[self.VRS_STARTS_FIELD].append(start)
                annotations[self.VRS_ENDS_FIELD].append(end)
                annotations[self.VRS_STATES_FIELD].append(alt)

    def _process_vcf_row(
        self,
        record: pysam.VariantRecord,
        vrs_data: dict,
        assembly: str,
        vrs_info_fields: list[str],
        incl_vrs_attrs: bool,
        incl_ref_allele: bool,
        create_pickle: bool,
        require_validation: bool,
    ) -> dict:
        """Compute VRS objects for a VCF row.

        Get VRS data for record's reference (if requested) and alt alleles. Return
        INFO field values to annotate VCF row with.

        # TODO update these
        """
        info_field_annotations = {field: [] for field in vrs_info_fields}

        # Get VRS data for reference allele
        gnomad_loc = f"{record.chrom}-{record.pos}"
        if incl_ref_allele:
            reference_allele = f"{gnomad_loc}-{record.ref}-{record.ref}"
            self._process_allele(
                reference_allele,
                vrs_data,
                info_field_annotations,
                assembly,
                create_pickle=create_pickle,
                incl_vrs_attrs=incl_vrs_attrs,
                require_validation=require_validation,
            )

        # Get VRS data for alts
        alts = record.alts or []
        alleles = [f"{gnomad_loc}-{record.ref}-{a}" for a in [*alts]]
        data_key = f"{record.chrom}\t{record.pos}\t{record.ref}\t{record.alts}"
        for allele in alleles:
            if "*" in allele:
                _logger.debug("Star allele found: %s", allele)
                for field in vrs_info_fields:
                    info_field_annotations[field].append("")
            else:
                self._process_allele(
                    allele,
                    vrs_data,
                    info_field_annotations,
                    assembly,
                    vrs_data_key=data_key,
                    create_pickle=create_pickle,
                    incl_vrs_attrs=incl_vrs_attrs,
                    require_validation=require_validation,
                )

        return info_field_annotations

    @use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.MISSING)
    def annotate(
        self,
        input_vcf_path: Path,
        output_vcf_path: Path | None = None,
        output_pkl_path: Path | None = None,
        incl_vrs_attrs: bool = False,
        incl_ref_allele: bool = True,
        assembly: str = "GRCh38",
        require_validation: bool = True,
    ) -> None:
        """Given a VCF, produce an output VCF annotated with VRS allele IDs, and/or
        a pickle file containing the full VRS objects.

        :param input_vcf_path: location of input VCF
        :param output_vcf_path: location at which to save output VCF (optional)
        :param output_pkl_path: location at which to save output PKL file (output)
        :param incl_vrs_attrs: whether ``VRS_Start``, ``VRS_End``, and ``VRS_State``
            attributes should be included in output VCF info field. These properties
            may be useful to retain outside of the VRS object for reasons like
            searchability. Does nothing if ``output_vcf_path`` left unset.
        :param incl_ref_allele: If true, perform VRS ID computation for REF allele and
            include the corresponding VRS object in any data dumps
        :param assembly: The assembly used in `vcf_in` data
        :param require_validation: If ``True``, validation checks (i.e., REF value
            matches expected REF at given location) must pass in order to return a VRS
            object for a record. If ``False`` then VRS object will be returned even if
            validation checks fail, although all instances of failed validation are
            logged as warnings regardless.
        """
        if not any((output_vcf_path, output_pkl_path)):
            msg = "Must provide one of: `output_vcf_path` or `output_pkl_path`"
            raise VCFAnnotatorError(msg)

        vcf = pysam.VariantFile(filename=str(input_vcf_path.absolute()))
        if output_vcf_path:
            self._update_vcf_header(vcf, incl_ref_allele, incl_vrs_attrs)

        vcf_out = (
            pysam.VariantFile(str(output_vcf_path.absolute()), "w", header=vcf.header)
            if output_vcf_path
            else None
        )
        create_pkl = bool(output_pkl_path)

        vrs_data = {}
        for record in vcf:
            if vcf_out:
                vrs_info_fields = [self.VRS_ALLELE_IDS_FIELD]
                if incl_vrs_attrs:
                    vrs_info_fields += [
                        self.VRS_STARTS_FIELD,
                        self.VRS_ENDS_FIELD,
                        self.VRS_STATES_FIELD,
                    ]
            else:
                # no info fields are necessary if we aren't producing an annotated VCF
                vrs_info_fields = []
            try:
                vrs_info_field_annotations = self._process_vcf_row(
                    record,
                    vrs_data,
                    assembly,
                    vrs_info_fields,
                    incl_vrs_attrs=incl_vrs_attrs,
                    incl_ref_allele=incl_ref_allele,
                    create_pickle=create_pkl,
                    require_validation=require_validation,
                )
            except Exception as ex:
                _logger.exception("VRS error on %s-%s", record.chrom, record.pos)
                err_msg = f"{ex}" or f"{type(ex)}"
                err_msg = err_msg.translate(self.VCF_ESCAPE_MAP)
                vrs_info_fields = [self.VRS_ERROR_FIELD]
                vrs_info_field_annotations = {self.VRS_ERROR_FIELD: [err_msg]}

            _logger.debug(
                "VCF record %s-%s generated vrs_field_data %s",
                record.chrom,
                record.pos,
                vrs_info_field_annotations,
            )
            if output_vcf_path and vcf_out:
                for k in vrs_info_fields:
                    record.info[k] = [
                        value or "." for value in vrs_info_field_annotations[k]
                    ]
                vcf_out.write(record)

        vcf.close()

        if vcf_out:
            vcf_out.close()
        if output_pkl_path:
            with output_pkl_path.open("wb") as wf:
                pickle.dump(vrs_data, wf)
