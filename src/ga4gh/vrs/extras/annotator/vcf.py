"""Annotate VCFs with VRS

$ vrs-annotate vcf input.vcf.gz --vcf_out output.vcf.gz --vrs_pickle_out vrs_objects.pkl

"""

import logging
import pickle
from pathlib import Path
from typing import ClassVar

import pysam
from pydantic import ValidationError

from ga4gh.core import VrsObjectIdentifierIs, use_ga4gh_compute_identifier_when
from ga4gh.vrs.dataproxy import DataProxyValidationError, _DataProxy
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

        :param data_proxy:
        """
        self.data_proxy = data_proxy
        self.tlr = AlleleTranslator(self.data_proxy)

    def _update_vcf_header(
        self,
        vcf: pysam.VariantFile,
        info_field_num: str,
        info_field_desc: str,
        incl_vrs_attrs: bool,
    ) -> None:
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

    @use_ga4gh_compute_identifier_when(VrsObjectIdentifierIs.MISSING)
    def annotate(
        self,
        input_vcf_path: Path,
        output_vcf_path: Path | None = None,
        output_pkl_path: Path | None = None,
        incl_vrs_attrs: bool = False,
        assembly: str = "GRCh38",
        compute_for_ref: bool = True,
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
        :param assembly: The assembly used in `vcf_in` data
        :param compute_for_ref: If true, compute VRS IDs for the reference allele
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
        info_field_num = "R" if compute_for_ref else "A"
        info_field_desc = "REF and ALT" if compute_for_ref else "ALT"
        self._update_vcf_header(vcf, info_field_num, info_field_desc, incl_vrs_attrs)

        vcf_out = (
            pysam.VariantFile(str(output_vcf_path.absolute()), "w", header=vcf.header)
            if output_vcf_path
            else None
        )

        vrs_data = {}
        for record in vcf:
            additional_info_fields = [self.VRS_ALLELE_IDS_FIELD]
            if incl_vrs_attrs:
                additional_info_fields += [
                    self.VRS_STARTS_FIELD,
                    self.VRS_ENDS_FIELD,
                    self.VRS_STATES_FIELD,
                ]
            try:
                vrs_field_data = self._get_vrs_data(
                    record,
                    vrs_data,
                    assembly,
                    additional_info_fields,
                    incl_vrs_attrs=incl_vrs_attrs,
                    output_pickle=output_pickle,
                    output_vcf=output_vcf,
                    compute_for_ref=compute_for_ref,
                    require_validation=require_validation,
                )
            except Exception as ex:
                _logger.exception("VRS error on %s-%s", record.chrom, record.pos)
                err_msg = f"{ex}" or f"{type(ex)}"
                err_msg = err_msg.translate(self.VCF_ESCAPE_MAP)
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
        vrs_data: dict,
        vrs_field_data: dict,
        assembly: str,
        vrs_data_key: str | None = None,
        output_pickle: bool = True,
        output_vcf: bool = False,
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
        :param output_pickle: If `True`, VRS pickle file will be output.
        :param output_vcf: If `True`, annotated VCF file will be output.
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

        if output_pickle and vrs_obj:
            key = vrs_data_key if vrs_data_key else vcf_coords
            vrs_data[key] = str(vrs_obj.model_dump(exclude_none=True))

        if output_vcf:
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
                    start = ""
                    end = ""
                    alt = ""

                vrs_field_data[self.VRS_STARTS_FIELD].append(start)
                vrs_field_data[self.VRS_ENDS_FIELD].append(end)
                vrs_field_data[self.VRS_STATES_FIELD].append(alt)

    def _get_vrs_data(
        self,
        record: pysam.VariantRecord,
        vrs_data: dict,
        assembly: str,
        additional_info_fields: list[str],
        incl_vrs_attrs: bool,
        output_pickle: bool,
        output_vcf: bool,
        compute_for_ref: bool,
        require_validation: bool,
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
        :param output_pickle: If `True`, VRS pickle file will be output.
        :param output_vcf: If `True`, annotated VCF file will be output.
        :param compute_for_ref: If true, compute VRS IDs for the reference allele
        :param require_validation: If `True` then validation checks must pass in
            order to return a VRS object. A `DataProxyValidationError` will be raised if
            validation checks fail. If `False` then VRS object will be returned even if
            validation checks fail. Defaults to `True`.
        :return: If `output_vcf = True`, a dictionary containing VRS Fields and list
            of associated values. If `output_vcf = False`, an empty dictionary will be
            returned.
        """
        vrs_field_data = (
            {field: [] for field in additional_info_fields} if output_vcf else {}
        )

        # Get VRS data for reference allele
        gnomad_loc = f"{record.chrom}-{record.pos}"
        if compute_for_ref:
            reference_allele = f"{gnomad_loc}-{record.ref}-{record.ref}"
            self._get_vrs_object(
                reference_allele,
                vrs_data,
                vrs_field_data,
                assembly,
                output_pickle=output_pickle,
                output_vcf=output_vcf,
                vrs_attributes=incl_vrs_attrs,
                require_validation=require_validation,
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
                    allele,
                    vrs_data,
                    vrs_field_data,
                    assembly,
                    vrs_data_key=data,
                    output_pickle=output_pickle,
                    output_vcf=output_vcf,
                    vrs_attributes=incl_vrs_attrs,
                    require_validation=require_validation,
                )

        return vrs_field_data
