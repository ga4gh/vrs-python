"""Translates various external formats into VRS models.

Input formats: VRS (serialized), hgvs, spdi, gnomad (vcf), beacon
Output formats: VRS (serialized), hgvs, spdi, gnomad (vcf)

"""
from typing import Tuple
from collections.abc import Mapping
import logging
import re

from bioutils.accessions import coerce_namespace
import hgvs.parser
import hgvs.location
import hgvs.posedit
import hgvs.edit
import hgvs.sequencevariant
import hgvs.dataproviders.uta

from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models, normalize as do_normalize
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
    gnomad_re = re.compile(r"(?P<chr>[^-]+)-(?P<pos>\d+)-(?P<ref>[ACGTN]+)-(?P<alt>[ACGTN]+|\*|\.)", re.IGNORECASE)
    hgvs_re = re.compile(r"[^:]+:[cgnpr]\.")
    spdi_re = re.compile(r"(?P<ac>[^:]+):(?P<pos>\d+):(?P<del_len_or_seq>\w*):(?P<ins_seq>\w*)")


    def __init__(self,  # pylint: disable=too-many-arguments
                 data_proxy,
                 default_assembly_name="GRCh38",
                 translate_sequence_identifiers=True,
                 normalize=True,
                 identify=True):
        self.default_assembly_name = default_assembly_name
        self.data_proxy = data_proxy
        self.translate_sequence_identifiers = translate_sequence_identifiers
        self.identify = identify
        self.normalize = normalize
        self.hgvs_tools = None


    def translate_from(self, var, fmt=None, **kwargs):
        """Translate variation `var` to VRS object


        If `fmt` is None, guess the appropriate format and return the variant.
        If `fmt` is specified, try only that format.

        kwargs:
            assembly_name (str): Assembly used for `var`. Defaults to the
                `default_assembly_name`. Only used for beacon and gnomad.
            require_validation (bool): If `True` then validation checks must pass in
                order to return a VRS object. A `ValidationError` will be raised if
                validation checks fail. If `False` then VRS object will be returned even
                if validation checks fail. Defaults to `True`.

        See also notes about `from_` and `to_` methods.
        """

        if fmt:
            t = self.from_translators[fmt]
            o = t(self, var, **kwargs)
            if o is None:
                raise ValueError(f"Unable to parse data as {fmt} variation")
            return o

        for _, t in self.from_translators.items():
            o = t(self, var, **kwargs)  # pylint: disable=too-many-function-args
            if o:
                return o

        formats = list(self.from_translators.keys())
        raise ValueError(f"Unable to parse data as {', '.join(formats)}")

    def translate_to(self, vo, fmt):
        """translate vrs object `vo` to named format `fmt`"""
        t = self.to_translators[fmt]
        return t(self, self.ensure_allele_is_latest_model(vo))

    def ensure_allele_is_latest_model(self, allele):
        """
        Change deprecated models:
        SequenceState -> LiteralSequenceExpression
        SimpleInterval -> SequenceInterval
        """
        if allele.state.type == "SequenceState":
            allele.state = models.LiteralSequenceExpression(
                sequence=allele.state.sequence,
            )
        if allele.location.interval.type == "SimpleInterval":
            allele.location.interval = models.SequenceInterval(
                start=models.Number(value=allele.location.interval.start),
                end=models.Number(value=allele.location.interval.end),
            )
        return allele


    ############################################################################
    ## INTERNAL

    def _from_beacon(self, beacon_expr, **kwargs):  # pylint: disable=too-many-locals
        """Parse beacon expression into VRS Allele

        kwargs:
            assembly_name (str): Assembly used for `beacon_expr`.

        #>>> a = tlr.from_beacon("19 : 44908822 C > T")
        #>>> a.as_dict()
        {'location': {'interval': {
           'end': {'value': 44908822, 'type': Number},
           'start': {'value': 44908821, 'type': Number},
           'type': 'SequenceInterval'},
          'sequence_id': 'GRCh38:19',
          'type': 'SequenceLocation'},
         'state': {'sequence': 'T', 'type': 'LiteralSequenceExpression'},
         'type': 'Allele'}

        """

        if not isinstance(beacon_expr, str):
            return None
        m = self.beacon_re.match(beacon_expr.replace(" ", ""))
        if not m:
            return None

        g = m.groupdict()
        assembly_name = kwargs.get("assembly_name", self.default_assembly_name)
        sequence_id = assembly_name + ":" + g["chr"]
        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        interval = models.SequenceInterval(start=models.Number(value=start),
                                           end=models.Number(value=end))
        location = models.SequenceLocation(sequence_id=sequence_id,
                                           interval=interval)
        state = models.LiteralSequenceExpression(sequence=ins_seq)
        allele = models.Allele(location=location, state=state)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_gnomad(self, gnomad_expr, **kwargs):  # pylint: disable=too-many-locals
        """Parse gnomAD-style VCF expression into VRS Allele

        :param str gnomad_expr: chr-pos-ref-alt

        kwargs:
            assembly_name (str): Assembly used for `gnomad_expr`.
            require_validation (bool): If `True` then validation checks must pass in
                order to return a VRS object. A `ValidationError` will be raised if
                validation checks fail. If `False` then VRS object will be returned even
                if validation checks fail. Defaults to `True`.

        #>>> a = tlr.from_gnomad("1-55516888-G-GA")
        #>>> a.as_dict()
        {'location': {'interval': {
           'end': {'value': 55516888, 'type': Number},
           'start': {'value': 55516887, 'type': Number},
           'type': 'SequenceInterval'},
          'sequence_id': 'GRCh38:1',
          'type': 'SequenceLocation'},
         'state': {'sequence': 'GA', 'type': 'LiteralSequenceExpression'},
         'type': 'Allele'}

        """
        if not isinstance(gnomad_expr, str):
            return None
        m = self.gnomad_re.match(gnomad_expr)
        if not m:
            return None

        g = m.groupdict()
        assembly_name = kwargs.get("assembly_name", self.default_assembly_name)
        sequence_id = assembly_name + ":" + g["chr"]
        start = int(g["pos"]) - 1
        ref = g["ref"].upper()
        alt = g["alt"].upper()
        end = start + len(ref)
        ins_seq = alt

        # validation checks
        valid_ref_seq, err_msg = self._is_valid_ref_seq(sequence_id, start, end, ref)
        if kwargs.get("require_validation", True) and not valid_ref_seq:
            raise ValidationError(err_msg)

        interval = models.SequenceInterval(start=models.Number(value=start),
                                           end=models.Number(value=end))
        location = models.SequenceLocation(sequence_id=sequence_id, interval=interval)
        sstate = models.LiteralSequenceExpression(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_hgvs(self, hgvs_expr, **kwargs):  # pylint: disable=unused-argument
        """parse hgvs into a VRS object (typically an Allele)

        #>>> a = tlr.from_hgvs("NM_012345.6:c.22A>T")
        #>>> a.as_dict()
        {
          'location': {
            'interval': {
               'end': {'value': 22, 'type': Number},
               'start': {'value': 21, 'type': Number},
               'type': 'SequenceInterval'},
            'sequence_id': 'refseq:NM_012345.6',
            'type': 'SequenceLocation'
          },
          'state': {'sequence': 'T', 'type': 'LiteralSequenceExpression'},
          'type': 'Allele'
        }

        """

        if not isinstance(hgvs_expr, str):
            return None
        if not self.hgvs_re.match(hgvs_expr):
            return None

        sv = self._hgvs_parser.parse_hgvs_variant(hgvs_expr)  # pylint: disable=no-member

        # prefix accession with namespace
        sequence_id = coerce_namespace(sv.ac)

        if isinstance(sv.posedit.pos, hgvs.location.BaseOffsetInterval):
            if sv.posedit.pos.start.is_intronic or sv.posedit.pos.end.is_intronic:
                raise ValueError("Intronic HGVS variants are not supported ({sv.posedit})")

        if sv.posedit.edit.type == "ins":
            interval = models.SequenceInterval(
                start=models.Number(value=sv.posedit.pos.start.base),
                end=models.Number(value=sv.posedit.pos.start.base))
            state = sv.posedit.edit.alt
        elif sv.posedit.edit.type in ("sub", "del", "delins", "identity"):
            interval = models.SequenceInterval(
                start=models.Number(value=sv.posedit.pos.start.base - 1),
                end=models.Number(value=sv.posedit.pos.end.base))
            if sv.posedit.edit.type == "identity":
                state = self.data_proxy.get_sequence(sv.ac,
                                                     sv.posedit.pos.start.base - 1,
                                                     sv.posedit.pos.end.base)
            else:
                state = sv.posedit.edit.alt or ""
        elif sv.posedit.edit.type == "dup":

            interval = models.SequenceInterval(
                start=models.Number(value=sv.posedit.pos.start.base - 1),
                end=models.Number(value=sv.posedit.pos.end.base))

            ref = self.data_proxy.get_sequence(sv.ac,
                                               sv.posedit.pos.start.base - 1,
                                               sv.posedit.pos.end.base)
            state = ref + ref
        else:
            raise ValueError(f"HGVS variant type {sv.posedit.edit.type} is unsupported")

        location = models.SequenceLocation(sequence_id=sequence_id, interval=interval)
        sstate = models.LiteralSequenceExpression(sequence=state)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_spdi(self, spdi_expr, **kwargs):  # pylint: disable=unused-argument
        """Parse SPDI expression in to a GA4GH Allele

        #>>> a = tlr.from_spdi("NM_012345.6:21:1:T")
        #>>> a.as_dict()
        {
          'location': {
            'interval': {
               'end': {'value': 22, 'type': Number},
               'start': {'value': 21, 'type': Number},
               'type': 'SequenceInterval'},
            'sequence_id': 'refseq:NM_012345.6',
            'type': 'SequenceLocation'
          },
          'state': {'sequence': 'T', 'type': 'LiteralSequenceExpression'},
          'type': 'Allele'
        }
        """

        if not isinstance(spdi_expr, str):
            return None
        m = self.spdi_re.match(spdi_expr)
        if not m:
            return None

        g = m.groupdict()
        sequence_id = coerce_namespace(g["ac"])
        start = int(g["pos"])
        try:
            del_len = int(g["del_len_or_seq"])
        except ValueError:
            del_len = len(g["del_len_or_seq"])
        end = start + del_len
        ins_seq = g["ins_seq"]

        interval = models.SequenceInterval(start=models.Number(value=start),
                                           end=models.Number(value=end))
        location = models.SequenceLocation(sequence_id=sequence_id, interval=interval)
        sstate = models.LiteralSequenceExpression(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_vrs(self, var, **kwargs):  # pylint: disable=unused-argument
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

    def _get_hgvs_tools(self):
        """ Only create UTA db connection if needed. There will be one connectionn per translator.
        """
        if self.hgvs_tools is None:
            uta_conn = hgvs.dataproviders.uta.connect()
            self.hgvs_tools = HgvsTools(uta_conn)
        return self.hgvs_tools

    def _to_hgvs(self, vo, namespace="refseq"):  # pylint: disable=too-many-locals
        """generates a *list* of HGVS expressions for VRS Allele.

        If `namespace` is not None, returns HGVS strings for the
        specified namespace.

        If `namespace` is None, returns HGVS strings for all alias
        translations.

        If no alias translations are available, an empty list is
        returned.

        If the VRS object cannot be expressed as HGVS, raises ValueError.

        """

        def ir_stype(a):
            if a.startswith("refseq:NM_"):
                return "n"
            if a.startswith("refseq:NP_"):
                return "p"
            if a.startswith("refseq:NG_"):
                return "g"
            if a.startswith("refseq:NC_"):
                return "g"
            if a.startswith("GRCh"):
                return "g"
            return None

        if not self.is_valid_allele(vo):
            raise ValueError("_to_hgvs requires a VRS Allele with SequenceLocation and LiteralSequenceExpression")

        sequence_id = str(vo.location.sequence_id)
        aliases = self.data_proxy.translate_sequence_identifier(sequence_id, namespace)

        # infer type of sequence based on accession
        # TODO: move to bioutils
        stypes = list(set(t for t in (ir_stype(a) for a in aliases) if t))
        if len(stypes) != 1:
            raise ValueError(f"Couldn't infer sequence type for {sequence_id} ({stypes})")
        stype = stypes[0]

        # build interval and edit depending on sequence type
        if stype == "p":  # pylint: disable=no-else-raise
            raise ValueError("Only nucleic acid variation is currently supported")
            # ival = hgvs.location.Interval(start=start, end=end)
            # edit = hgvs.edit.AARefAlt(ref=None, alt=vo.state.sequence)
        else:
            start, end = vo.location.interval.start.value, vo.location.interval.end.value
            # ib: 0 1 2 3 4 5
            #  h:  1 2 3 4 5
            if start == end:    # insert: hgvs uses *exclusive coords*
                ref = None
                end += 1
            else:               # else: hgvs uses *inclusive coords*
                ref = self.data_proxy.get_sequence(sequence_id, start, end)
                start += 1
            ival = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(base=start),
                end=hgvs.location.SimplePosition(base=end))
            alt = str(vo.state.sequence) or None  # "" => None
            edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)

        posedit = hgvs.posedit.PosEdit(pos=ival, edit=edit)
        var = hgvs.sequencevariant.SequenceVariant(
            ac=None,
            type=stype,
            posedit=posedit)

        hgvs_exprs = []
        for alias in aliases:
            ns, a = alias.split(":")
            # skip GRCh accessions unless specifically requested
            # because they are ambiguous without their namespace,
            # which can't be included in HGVS expressions
            # TODO: use default_assembly_name here
            if ns.startswith("GRC") and namespace is None:
                continue

            if not any(a.startswith(pfx) for pfx in ("NM", "NP", "NC", "NG")):
                continue

            var.ac = a

            try:
                if not namespace.startswith("GRC"):
                    # if the namespace is GRC, can't normalize, since hgvs can't deal with it
                    hgvs_tools = self._get_hgvs_tools()
                    parsed = hgvs_tools.parse(str(var))
                    var = hgvs_tools.normalize(parsed)

                hgvs_exprs += [str(var)]
            except hgvs.exceptions.HGVSDataNotAvailableError:
                _logger.warning("No data found for accession %s", a)

        return list(set(hgvs_exprs))


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
        if not self.is_valid_allele(vo):
            raise ValueError("_to_spdi requires a VRS Allele with SequenceLocation and LiteralSequenceExpression")

        sequence_id = str(vo.location.sequence_id)
        aliases = self.data_proxy.translate_sequence_identifier(sequence_id, namespace)
        aliases = [a.split(":")[1] for a in aliases]
        start, end = vo.location.interval.start.value, vo.location.interval.end.value
        spdi_tail = f":{start}:{end-start}:{vo.state.sequence}"
        spdis = [a + spdi_tail for a in aliases]
        return spdis


    def _is_valid_ref_seq(self, sequence_id: str, start_pos: int, end_pos: int,
                          ref: str) -> Tuple[bool, str]:
        """Return wether or not the expected reference sequence matches the actual reference sequence

        :param str sequence_id: Sequence ID to use
        :param int start_pos: Start pos (inter-residue) on the sequence_id
        :param int end_pos: End pos (inter-residue) on the sequence_id
        :param str ref: The expected reference sequence on the sequence_id given the
            start_pos and end_pos
        :return: Tuple containing whether or not actual reference sequence matches
            the expected reference sequence and error message if mismatch
        """
        actual_ref = self.data_proxy.get_sequence(sequence_id, start_pos, end_pos)
        is_valid = actual_ref == ref
        err_msg = ""
        if not is_valid:
            err_msg = f"Expected reference sequence {ref} on {sequence_id} at positions "\
                      f"({start_pos}, {end_pos}) but found {actual_ref}"
            _logger.warning(err_msg)
        return is_valid, err_msg


    def is_valid_allele(self, vo):
        """Ensure that `vo` is a valid VRS Allele with SequenceLocation and
        LiteralSequenceExpression
        """
        return (vo.type == "Allele"
                and vo.location.type == "SequenceLocation"
                and vo.state.type == "LiteralSequenceExpression")


    @lazy_property
    def _hgvs_parser(self):
        """instantiates and returns an hgvs parser instance"""
        _logger.info("Creating  parser")
        return hgvs.parser.Parser()


    def _post_process_imported_allele(self, allele):
        """Provide common post-processing for imported Alleles IN-PLACE.

        """

        if self.translate_sequence_identifiers:
            seq_id = self.data_proxy.translate_sequence_identifier(allele.location.sequence_id._value, "ga4gh")[0]
            allele.location.sequence_id = seq_id

        if self.normalize:
            allele = do_normalize(allele, self.data_proxy)

        if self.identify:
            allele._id = ga4gh_identify(allele)
            allele.location._id = ga4gh_identify(allele.location)

        return allele


    def _seq_id_mapper(self, ir):
        if self.translate_sequence_identifiers:
            return self.data_proxy.translate_sequence_identifier(ir, "ga4gh")[0]
        return ir


    from_translators = {
        "beacon": _from_beacon,
        "gnomad": _from_gnomad,
        "hgvs": _from_hgvs,
        "spdi": _from_spdi,
        "vrs": _from_vrs,
    }

    to_translators = {
        "hgvs": _to_hgvs,
        "spdi": _to_spdi,
        #"gnomad": to_gnomad,
    }




if __name__ == "__main__":
    # pylint: disable=ungrouped-imports

    import coloredlogs
    coloredlogs.install(level="INFO")

    from ga4gh.vrs.dataproxy import create_dataproxy
    dp = create_dataproxy("seqrepo+file:///usr/local/share/seqrepo/latest")
    tlr = Translator(data_proxy=dp)

    expressions = [
        "bogus",
        "1-55516888-G-GA",
        "19 : 44908822 C > T",
        "NC_000019.10:g.44908822C>T",
        "NM_000551.3:21:1:T", {
            "location": {
                "interval": {
                    "end": {"value": 22, "type": "Number"},
                    "start": {"value": 21, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        }, {
           "end": {"value": 22, "type": "Number"},
            "start": {"value": 21, "type": "Number"},
            "type": "SequenceInterval"
        }
    ]
    from_formats = ["hgvs", "gnomad", "beacon", "spdi", "vrs", None]

    for e in expressions:
        print(f"* {e}")
        for f in from_formats:
            try:
                vrs_obj = tlr.translate_from(e, f)
                r = vrs_obj.type
            except ValueError:
                r = "-"
            except Exception as ex:  # pylint: disable=broad-except
                r = ex.__class__.__name__
            print(f"  {f}: {r}")
