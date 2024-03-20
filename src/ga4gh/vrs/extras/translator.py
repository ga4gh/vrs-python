"""Translates various external formats into VRS models.

Input formats: VRS (serialized), hgvs, spdi, gnomad (vcf), beacon
Output formats: VRS (serialized), hgvs, spdi, gnomad (vcf)

"""

from collections.abc import Mapping
from typing import Optional, Union, Tuple
from vrs import _DataProxy
import logging
import re

from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models, normalize
from ga4gh.vrs.extras.decorators import lazy_property  # this should be relocated
from ga4gh.vrs.utils.hgvs_tools import HgvsTools

_logger = logging.getLogger(__name__)


class ValidationError(Exception):
    """Class for validation errors during translation"""


class Translator:
    """Translates various variation formats to and from GA4GH VRS models

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
    hgvs_re = re.compile(r"[^:]+:[cgmnpr]\.")
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
        self.hgvs_tools = None
        self.from_translators = {}
        self.to_translators = {}
    
    # def _is_intronic(self, sv: hgvs.sequencevariant.SequenceVariant):
    #     """Check if HGVS sequence variant is intronic"""
    #     if isinstance(sv.posedit.pos, hgvs.location.BaseOffsetInterval):
    #         if sv.posedit.pos.start.is_intronic or sv.posedit.pos.end.is_intronic:
    #             return True
    #     return False

    # def _get_hgvs_refget_ac(self, sv: hgvs.sequencevariant.SequenceVariant):
    #     """Get refget accession for hgvs sequence variant"""
    #     # prefix accession with namespace
    #     sequence = coerce_namespace(sv.ac)
    #     return self._get_refget_accession(sequence)

    @staticmethod
    def _ir_stype(a):
        """
        The purpose of this function is to provide a convenient way to extract the sequence type from an accession by matching its prefix to a known set of prefixes.

        Args:
        a (str): The accession string.

        Returns:
        str or None: The sequence type associated with the accession string, or None if no matching prefix is found.
        """

        prefix_dict = {
            "refseq:NM_": "n",
            "refseq:NC_012920": "m",
            "refseq:NG_": "g",
            "refseq:NC_00": "g",
            "refseq:NW_": "g",
            "refseq:NR_": "n",
            "refseq:NP_": "p",
            "refseq:XM_": "n",
            "refseq:XR_": "n",
            "refseq:XP_": "p",
            "GRCh": "g",
        }

        for prefix, stype in prefix_dict.items():
            if a.startswith(prefix):
                return stype
        return None

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
                    order to return a VRS object. A `ValidationError` will be raised if
                    validation checks fail. If `False` then VRS object will be returned
                    even if validation checks fail. Defaults to `True`.
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

    def _get_hgvs_tools(self):
        """ Only create UTA db connection if needed. There will be one connection per
        translator.
        """
        if self.hgvs_tools is None:
            self.hgvs_tools = HgvsTools(self.data_proxy)
        return self.hgvs_tools

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

class AlleleTranslator(Translator):
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

    def _create_allele(self, refget_accession, start, end, ins_seq, **kwargs):
        """
        Create an allele object with the given parameters.

        Args:
            refget_accession (str): The accession ID of the reference genome.
            start (int): The start position of the allele.
            end (int): The end position of the allele.
            ins_seq (str): The inserted sequence for the allele.
            **kwargs: Additional keyword arguments.

        Returns:
            models.Allele: The created allele object.
        """
        seq_ref = models.SequenceReference(refgetAccession=refget_accession)
        location = models.SequenceLocation(sequenceReference=seq_ref, start=start, end=end)
        state = models.LiteralSequenceExpression(sequence=ins_seq)
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
        refget_accession = self.data_proxy.get_refget_accession(sequence)
        if not refget_accession:
            return None

        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        allele = self._create_allele(refget_accession, start, end, ins_seq, **kwargs)

        return allele


    def _from_gnomad(self, gnomad_expr, **kwargs):
        """Parse gnomAD-style VCF expression into VRS Allele

        kwargs:
            assembly_name (str): Assembly used for `gnomad_expr`.
            require_validation (bool): If `True` then validation checks must pass in
                order to return a VRS object. A `ValidationError` will be raised if
                validation checks fail. If `False` then VRS object will be returned even
                if validation checks fail. Defaults to `True`.
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
        refget_accession = self.data_proxy.get_refget_accession(sequence)
        # refget_accession = self._get_refget_accession(sequence)
        if not refget_accession:
            return None

        start = int(g["pos"]) - 1
        ref = g["ref"].upper()
        alt = g["alt"].upper()
        end = start + len(ref)
        ins_seq = alt

        # TODO: ask other devs if this should be down on all _from_...  methods?
        if kwargs.get("require_validation", True):
            # validation check for matching reference sequence bases
            valid_ref_seq, err_msg = self.data_proxy.is_valid_ref_seq(sequence, start, end, ref)
            if not valid_ref_seq:
                raise ValidationError(err_msg)

        allele = self._create_allele(refget_accession, start, end, ins_seq, **kwargs)
        return allele

    def _from_hgvs(self, hgvs_expr: str, **kwargs):
        """parse hgvs into a VRS Allele Object

        kwargs:
            rle_seq_limit Optional(int): If RLE is set as the new state after
                normalization, this sets the limit for the length of the `sequence`.
                To exclude `sequence` from the response, set to 0.
                For no limit, set to `None`.
                Defaults value set in instance variable, `rle_seq_limit`.
            do_normalize (bool): `True` if fully justified normalization should be
                performed. `False` otherwise. Defaults to `True`

        #>>> a = tlr.from_hgvs("NC_000007.14:g.55181320A>T")
        #>>> a.model_dump()
        {
          'location': {
            'end': 55181320,
            'start': 55181319,
            'sequenceReference': {
              'type': 'SequenceReference',
              'refgetAccession': 'SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul'
            },
            'type': 'SequenceLocation'
          },
          'state': {
            'sequence': 'T',
            'type': 'LiteralSequenceExpression'
          },
          'type': 'Allele'
        }

        """
        sv = self._get_hgvs_tools().parse(hgvs_expr)
        if not sv:
            return None
        
        if self._get_hgvs_tools().is_intronic(sv):
            raise ValueError("Intronic HGVS variants are not supported")

        refget_accession = self.data_proxy.derive_refget_accession(sv.ac)
        if not refget_accession:
            return None
        
        # translate coding coordinates to positional coordinates, if necessary
        if sv.type == "c":
            sv = self._get_hgvs_tools().c_to_n(sv)

        # TODO: shouldn't we be validating the "ref" sequence here? if the kwargs "require_validation" value is set to True

        (start,end,state) = self._get_hgvs_tools().get_position_and_state(sv)

        allele = self._create_allele(refget_accession, start, end, state, **kwargs)
        return allele

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
        refget_accession = self.data_proxy.get_refget_accession(g["ac"])
        # # TODO #1: move the coerce_namespace into the _get_refget_accession function to consolidate the logic
        # # TODO #2: move _get_refget_accession into data_proxy
        # sequence = coerce_namespace(g["ac"])
        # refget_accession = self._get_refget_accession(sequence)
        if not refget_accession:
            return None

        start = int(g["pos"])
        try:
            del_len = int(g["del_len_or_seq"])
        except ValueError:
            del_len = len(g["del_len_or_seq"])
        end = start + del_len
        ins_seq = g["ins_seq"]

        allele = self._create_allele(refget_accession, start, end, ins_seq, **kwargs)
        return allele

    def _to_hgvs(self, vo, namespace="refseq"):
        """generates a *list* of HGVS expressions for VRS Allele.

        If `namespace` is not None, returns HGVS strings for the
        specified namespace.

        If `namespace` is None, returns HGVS strings for all alias
        translations.

        If no alias translations are available, an empty list is
        returned.

        If the VRS object cannot be expressed as HGVS, raises ValueError.

        This method assumes that IRIs are dereferenced, providing a `SequenceReference`
        as the `vo.location.sequenceReference`. If a `SequenceReference` is not
        provided, raises TypeError
        """


        # # TODO: this method should extract the VRS elements needed to build the HGVS expression 
        # # and pass them to a method in the hgvs_tools class to build the final result
        return None
    
        # if not isinstance(vo.location.sequenceReference, models.SequenceReference):
        #     raise TypeError(
        #         "`vo.location.sequenceReference` expects a `SequenceReference`"
        #     )

        # sequence = f"ga4gh:{export_sequencelocation_sequence_id(vo.location.sequenceReference)}"
        # aliases = self.data_proxy.translate_sequence_identifier(sequence, namespace)

        # # infer type of sequence based on accession
        # # TODO: move to bioutils
        # stypes = list(set(t for t in (self._ir_stype(a) for a in aliases) if t))
        # if len(stypes) != 1:
        #     # TODO: consider returning all stypes if more than one alias type is returned
        #     raise ValueError(f"Couldn't infer sequence type for {sequence} ({stypes})")
        # stype = stypes[0]

        # # build interval and edit depending on sequence type
        # if stype == "p":
        #     raise ValueError("Only nucleic acid variation is currently supported")
        #     # ival = hgvs.location.Interval(start=start, end=end)
        #     # edit = hgvs.edit.AARefAlt(ref=None, alt=vo.state.sequence)
        # else:                   # pylint: disable=no-else-raise
        #     start, end = vo.location.start, vo.location.end
        #     # ib: 0 1 2 3 4 5
        #     #  h:  1 2 3 4 5
        #     if start == end:    # insert: hgvs uses *exclusive coords*
        #         ref = None
        #         end += 1
        #     else:               # else: hgvs uses *inclusive coords*
        #         ref = self.data_proxy.get_sequence(sequence, start, end)
        #         start += 1

        #     ival = hgvs.location.Interval(
        #         start=hgvs.location.SimplePosition(start),
        #         end=hgvs.location.SimplePosition(end)
        #     )
        #     alt = str(vo.state.sequence.root) or None  # "" => None
        #     edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)

        # posedit = hgvs.posedit.PosEdit(pos=ival, edit=edit)
        # var = hgvs.sequencevariant.SequenceVariant(
        #     ac=None,
        #     type=stype,
        #     posedit=posedit)

        # hgvs_exprs = []
        # for alias in aliases:
        #     ns, a = alias.split(":")
        #     # skip GRCh accessions unless specifically requested
        #     # because they are ambiguous without their namespace,
        #     # which can't be included in HGVS expressions
        #     # TODO: use default_assembly_name here
        #     if ns.startswith("GRC") and namespace is None:
        #         continue

        #     if not (any(a.startswith(pfx) for pfx in ("NM", "NP", "NC", "NG", "NR", "NW", "XM", "XR", "XP"))):
        #         continue

        #     var.ac = a

        #     try:
        #         if not namespace.startswith("GRC"):
        #             # if the namespace is GRC, can't normalize, since hgvs can't deal with it
        #             hgvs_tools = self._get_hgvs_tools()
        #             parsed = hgvs_tools.parse(str(var))
        #             var = hgvs_tools.normalize(parsed)

        #         hgvs_exprs += [str(var)]
                
        #     except hgvs.exceptions.HGVSDataNotAvailableError:
        #         _logger.warning(f"No data found for accession {a}")

        # return list(set(hgvs_exprs))

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
        sequence = f"ga4gh:{export_sequencelocation_sequence_id(vo.location.sequenceReference)}"
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


def export_sequencelocation_sequence_id(
    location_sequence_reference: Union[models.IRI, models.SequenceReference]
):
    if isinstance(location_sequence_reference, models.IRI):
        return location_sequence_reference.root
    elif isinstance(location_sequence_reference, models.SequenceReference):
        return location_sequence_reference.refgetAccession


class CnvTranslator(Translator):
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
        sv = self._get_hgvs_tools().parse(hgvs_dup_del_expr)
        if not sv:
            return None
        
        if self._get_hgvs_tools().get_edit_type(sv) not in {"del", "dup"}:
            raise ValueError("Must provide a 'del' or 'dup'")
       
        if self._get_hgvs_tools().is_intronic(sv):
            raise ValueError("Intronic HGVS variants are not supported")

        refget_accession = self.data_proxy.derive (sv)
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

if __name__ == "__main__":
    # pylint: disable=ungrouped-imports

    import coloredlogs
    coloredlogs.install(level="INFO")

    from ga4gh.vrs.dataproxy import create_dataproxy
from dataproxy import _DataProxy
    # dp = create_dataproxy("seqrepo+file:///usr/local/share/seqrepo/latest")
    dp = create_dataproxy("seqrepo + http://localhost:5000/seqrepo")
    tlr = Translator(data_proxy=dp)

    expressions = [
        "bogus",
        "1-55516888-G-GA",
        "19 : 44908822 C > T",
        "NC_000019.10:g.44908822C>T",
        "NM_000551.3:21:1:T", {
            "location": {
                "end": 22,
                "start": 21,
                "sequenceReference": {
                    "refgetAccession": "SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_",
                    "type": "SequenceReference"
                },
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        }, {
            "end": 22,
            "start": 21,
        }
    ]
    formats = ["hgvs", "gnomad", "beacon", "spdi", "vrs", None]

    for e in expressions:
        print(f"* {e}")
        for f in formats:
            try:
                o = tlr.translate_from(e, f)
                r = o.type
            except ValueError:
                r = "-"
            except Exception as ex:
                r = ex.__class__.__name__
            print(f"  {f}: {r}")
