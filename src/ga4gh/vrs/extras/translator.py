"""Translates various external formats into VRS models.

Input formats: VRS (serialized), hgvs, spdi, gnomad (vcf), beacon
Output formats: VRS (serialized), hgvs, spdi, gnomad (vcf)

"""

from abc import ABC
from collections.abc import Mapping
from typing import Optional, Union
from ga4gh.vrs.dataproxy import create_dataproxy, _DataProxy
from ga4gh.vrs.extras.decorators import lazy_property
import logging
import re

from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models, normalize
from ga4gh.vrs.utils.hgvs_tools import HgvsTools

_logger = logging.getLogger(__name__)

class _Translator(ABC):
    """abstract class / interface for VRS to/from translation needs

     Translates various variation formats to and from GA4GH VRS models

    All `from_` methods follow this pattern:
    * If the argument does not appear to be an appropriate type, None is returned
    * Otherwise, the argument is expected to be of the correct type.  If an error occurs during processing,
      an exception is raised.
    * Otherwise, the VRS object is returned

    """

    beacon_re = re.compile(r"(?P<chr>[^-]+)\s*:\s*(?P<pos>\d+)\s*(?P<ref>\w+)\s*>\s*(?P<alt>\w+)")
    gnomad_re = re.compile(
        r"(?P<chr>[^-]+)-(?P<pos>\d+)-(?P<ref>[ACGTURYKMSWBDHVN]+)-(?P<alt>[ACGTURYKMSWBDHVN]+)",
        re.IGNORECASE
    )
    spdi_re = re.compile(r"(?P<ac>[^:]+):(?P<pos>\d+):(?P<del_len_or_seq>\w*):(?P<ins_seq>\w*)")

    def __init__(
        self,
        data_proxy: _DataProxy,
        default_assembly_name="GRCh38",
        identify=True,
        rle_seq_limit: Optional[int] = 50
    ):
        self.default_assembly_name = default_assembly_name
        self.data_proxy = data_proxy
        self.identify = identify
        self.rle_seq_limit = rle_seq_limit
        self.from_translators = {}
        self.to_translators = {}
        return

    def translate_from(self, var, fmt=None, **kwargs):
        """Translate variation `var` to VRS object

        If `fmt` is None, guess the appropriate format and return the variant.
        If `fmt` is specified, try only that format.

        See also notes about `from_` and `to_` methods.

        kwargs:
            For CnvTranslator
                copies(int): The number of copies to use. If provided will return a
                    CopyNumberCount
                copy_change(models.CopyChange): Copy change. If not provided, default is
                    efo:0030067 for deletions and efo:0030070 for duplications
            For AlleleTranslator
                assembly_name (str): Assembly used for `var`. Defaults to the
                    `default_assembly_name`. Only used for beacon and gnomad.
                require_validation (bool): If `True` then validation checks must pass in
                    order to return a VRS object. A `DataProxyValidationError` will be
                    raised if validation checks fail. If `False` then VRS object will be
                    returned even if validation checks fail. Defaults to `True`.
                rle_seq_limit Optional(int): If RLE is set as the new state after
                    normalization, this sets the limit for the length of the `sequence`.
                    To exclude `sequence` from the response, set to 0.
                    For no limit, set to `None`.
                    Defaults value set in instance variable, `rle_seq_limit`.
                do_normalize (bool): `True` if fully justified normalization should be
                    performed. `False` otherwise. Defaults to `True`
        """
        if fmt:
            try:
                t = self.from_translators[fmt]
            except KeyError:
                raise NotImplementedError(f"{fmt} is not supported")
            else:
                o = t(var, **kwargs)
                if o is None:
                    raise ValueError(f"Unable to parse data as {fmt} variation")
                return o

        for _, t in self.from_translators.items():
            o = t(var, **kwargs)
            if o:
                return o

        formats = list(self.from_translators.keys())
        raise ValueError(f"Unable to parse data as {', '.join(formats)}")

    def translate_to(self, vo, fmt):
        """translate vrs object `vo` to named format `fmt`"""
        t = self.to_translators[fmt]
        return t(vo)

    ############################################################################
    # INTERNAL

    @lazy_property
    def hgvs_tools(self):
        """instantiates and returns an HgvsTools instance"""
        return HgvsTools(self.data_proxy)

    def _from_vrs(self, var):
        """convert from dict representation of VRS JSON to VRS object"""
        if not isinstance(var, Mapping):
            return None
        if "type" not in var:
            return None
        try:
            model = models[var["type"]]
        except KeyError:
            return None
        return model(**var)

class AlleleTranslator(_Translator):
    """Class for translating formats to and from VRS Alleles"""

    def __init__(
        self, data_proxy, default_assembly_name="GRCh38", identify=True
    ):
        """Initialize AlleleTranslator class"""
        super().__init__(data_proxy, default_assembly_name, identify)

        self.from_translators = {
            "beacon": self._from_beacon,
            "gnomad": self._from_gnomad,
            "hgvs": self._from_hgvs,
            "spdi": self._from_spdi,
            "vrs": self._from_vrs,
        }

        self.to_translators = {
            "hgvs": self._to_hgvs,
            "spdi": self._to_spdi,
        }

    def _create_allele(self, values: dict, **kwargs):
        """
        Create an allele object with the given parameters.

        Args:
            values (dict): The values to use for creating the allele object.
                'refget_accession' (str): The accession ID of the reference genome.
                'start' (int): The start position of the allele.
                'end' (int): The end position of the allele.
                literal_sequence' (str): The literal sequence for the allele.
            **kwargs: Additional keyword arguments.

        Returns:
            models.Allele: The created allele object.
        """
        seq_ref = models.SequenceReference(refgetAccession=values["refget_accession"])
        location = models.SequenceLocation(sequenceReference=seq_ref, start=values["start"], end=values["end"])
        state = models.LiteralSequenceExpression(sequence=values["literal_sequence"])
        allele = models.Allele(location=location, state=state)
        allele = self._post_process_imported_allele(allele, **kwargs)
        return allele

    def _from_beacon(self, beacon_expr, **kwargs):
        """Parse beacon expression into VRS Allele

        kwargs:
            assembly_name (str): Assembly used for `beacon_expr`.
            rle_seq_limit Optional(int): If RLE is set as the new state after
                normalization, this sets the limit for the length of the `sequence`.
                To exclude `sequence` from the response, set to 0.
                For no limit, set to `None`.
                Defaults value set in instance variable, `rle_seq_limit`.
            do_normalize (bool): `True` if fully justified normalization should be
                performed. `False` otherwise. Defaults to `True`

        #>>> a = tlr.from_beacon("19 : 44908822 C > T")
        #>>> a.model_dump()
        {
          'location': {
            'end': 44908822,
            'start': 44908821,
            'sequenceReference': {
              'type': 'SequenceReference',
              'refgetAccession': 'SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl'
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

        if not isinstance(beacon_expr, str):
            return None

        m = self.beacon_re.match(beacon_expr.replace(" ", ""))
        if not m:
            return None

        g = m.groupdict()
        assembly_name = kwargs.get("assembly_name", self.default_assembly_name)
        sequence = assembly_name + ":" + g["chr"]
        refget_accession = self.data_proxy.derive_refget_accession(sequence)
        if not refget_accession:
            return None

        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        values = {"refget_accession": refget_accession, "start": start, "end": end, "literal_sequence": ins_seq}
        allele = self._create_allele(values, **kwargs)

        return allele

    def _from_gnomad(self, gnomad_expr, **kwargs):
        """Parse gnomAD-style VCF expression into VRS Allele

        kwargs:
            assembly_name (str): Assembly used for `gnomad_expr`.
            require_validation (bool): If `True` then validation checks must pass in
                order to return a VRS object. A `DataProxyValidationError` will be
                raised if validation checks fail. If `False` then VRS object will be
                returned even if validation checks fail. Defaults to `True`.
            rle_seq_limit Optional(int): If RLE is set as the new state after
                normalization, this sets the limit for the length of the `sequence`.
                To exclude `sequence` from the response, set to 0.
                For no limit, set to `None`.
                Defaults value set in instance variable, `rle_seq_limit`.
            do_normalize (bool): `True` if fully justified normalization should be
                performed. `False` otherwise. Defaults to `True`

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
        m = self.gnomad_re.match(gnomad_expr)
        if not m:
            return None

        g = m.groupdict()
        assembly_name = kwargs.get("assembly_name", self.default_assembly_name)
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
            require_validation=kwargs.get("require_validation", True)
        )

        values = {"refget_accession": refget_accession, "start": start, "end": end, "literal_sequence": ins_seq}
        allele = self._create_allele(values, **kwargs)
        return allele

    def _from_hgvs(self, hgvs_expr: str, **kwargs):
        allele_values = self.hgvs_tools.extract_allele_values(hgvs_expr)
        if allele_values:
            return self._create_allele(allele_values, **kwargs)
        else:
            return None

    def _from_spdi(self, spdi_expr, **kwargs):
        """Parse SPDI expression in to a GA4GH Allele

        kwargs:
            rle_seq_limit Optional(int): If RLE is set as the new state after
                normalization, this sets the limit for the length of the `sequence`.
                To exclude `sequence` from the response, set to 0.
                For no limit, set to `None`.
                Defaults value set in instance variable, `rle_seq_limit`.
            do_normalize (bool): `True` if fully justified normalization should be
                performed. `False` otherwise. Defaults to `True`

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

        m = self.spdi_re.match(spdi_expr)
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

        values = {"refget_accession": refget_accession, "start": start, "end": end, "literal_sequence": ins_seq}

        allele = self._create_allele(values, **kwargs)
        return allele

    def _to_hgvs(self, vo, namespace="refseq"):
        return self.hgvs_tools.from_allele(vo, namespace)

    def _to_spdi(self, vo, namespace="refseq"):
        """generates a *list* of SPDI expressions for VRS Allele.

        If `namespace` is not None, returns SPDI strings for the
        specified namespace.

        If `namespace` is None, returns SPDI strings for all alias
        translations.

        If no alias translations are available, an empty list is
        returned.

        If the VRS object cannot be expressed as SPDI, raises ValueError.

        SPDI and VRS use identical normalization. The incoming Allele
        is expected to be normalized per VRS spec.

        """
        sequence = f"ga4gh:{vo.location.get_refget_accession()}"
        aliases = self.data_proxy.translate_sequence_identifier(sequence, namespace)
        aliases = [a.split(":")[1] for a in aliases]
        start, end = vo.location.start, vo.location.end
        spdi_tail = f":{start}:{end-start}:{vo.state.sequence.root}"
        spdis = [a + spdi_tail for a in aliases]
        return spdis

    def _post_process_imported_allele(self, allele, **kwargs):
        """Provide common post-processing for imported Alleles IN-PLACE.

        :param allele: VRS Allele object

        kwargs:
            rle_seq_limit Optional(int): If RLE is set as the new state after
                normalization, this sets the limit for the length of the `sequence`.
                To exclude `sequence` from the response, set to 0.
                For no limit, set to `None`.
            do_normalize (bool): `True` if fully justified normalization should be
                performed. `False` otherwise. Defaults to `True`
        """
        if kwargs.get("do_normalize", True):
            allele = normalize(
                allele,
                self.data_proxy,
                rle_seq_limit=kwargs.get("rle_seq_limit", self.rle_seq_limit)
            )

        if self.identify:
            allele.id = ga4gh_identify(allele)
            allele.location.id = ga4gh_identify(allele.location)

        return allele


class CnvTranslator(_Translator):
    """Class for translating formats from format to VRS Copy Number"""

    def __init__(
        self, data_proxy, default_assembly_name="GRCh38", identify=True
    ):
        """Initialize CnvTranslator class"""
        super().__init__(data_proxy, default_assembly_name, identify)
        self.from_translators = {
            "hgvs": self._from_hgvs,
        }

    def _from_hgvs(self, hgvs_dup_del_expr: str, **kwargs):
        """parse hgvs into a VRS CNV Object

        kwargs:
            copies: The number of copies to use. If provided will return a
                CopyNumberCount
            copy_change: Copy change. If not provided, default is efo:0030067 for
                deletions and efo:0030070 for duplications
        """
        # sv = self._get_parsed_hgvs(hgvs_dup_del_expr)
        sv = self.hgvs_tools.parse(hgvs_dup_del_expr)
        if not sv:
            return None

        sv_type = self.hgvs_tools.get_edit_type(sv)
        if sv_type not in {"del", "dup"}:
            raise ValueError("Must provide a 'del' or 'dup'")

        if self.hgvs_tools.is_intronic(sv):
            raise ValueError("Intronic HGVS variants are not supported")

        refget_accession = self.data_proxy.derive_refget_accession(sv.ac)
        if not refget_accession:
            return None

        location = models.SequenceLocation(
            sequenceReference=models.SequenceReference(refgetAccession=refget_accession),
            start=sv.posedit.pos.start.base - 1,
            end=sv.posedit.pos.end.base
        )

        copies = kwargs.get("copies")
        if copies:
            cnv = models.CopyNumberCount(location=location, copies=copies)
        else:
            copy_change = kwargs.get("copy_change")
            if not copy_change:
                copy_change = models.CopyChange.EFO_0030067 if sv_type == "del" else models.CopyChange.EFO_0030070
            cnv = models.CopyNumberChange(location=location, copyChange=copy_change)

        cnv =self._post_process_imported_cnv(cnv)
        return cnv

    def _post_process_imported_cnv(self, copy_number):
        """Provide common post-processing for imported Copy Numbers IN-PLACE."""
        if self.identify:
            copy_number.id = ga4gh_identify(copy_number)
            copy_number.location.id = ga4gh_identify(copy_number.location)

        return copy_number
