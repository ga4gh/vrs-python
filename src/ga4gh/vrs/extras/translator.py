"""Translates various external formats into VRS models.

Input formats: VRS (serialized), hgvs, spdi, gnomad (vcf), beacon
Output formats: VRS (serialized), hgvs, spdi, gnomad (vcf)

"""

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
from ga4gh.vrs import models, normalize
from ga4gh.vrs.extras.decorators import lazy_property  # this should be relocated
from ga4gh.vrs.utils.hgvs_tools import HgvsTools

_logger = logging.getLogger(__name__)


class Translator:
    """Translates various variation formats to and from GA4GH VRS models

    All `from_` methods follow this pattern:
    * If the argument does not appear to be an appropriate type, None is returned
    * Otherwise, the argument is expected to be of the correct type.  If an error occurs during processing,
      an exception is raised.
    * Otherwise, the VRS object is returned

    """

    beacon_re = re.compile(r"(?P<chr>[^-]+)\s*:\s*(?P<pos>\d+)\s*(?P<ref>\w+)\s*>\s*(?P<alt>\w+)")
    gnomad_re = re.compile(r"(?P<chr>[^-]+)-(?P<pos>\d+)-(?P<ref>\w+)-(?P<alt>\w+)")
    hgvs_re = re.compile(r"[^:]+:[cgnpr]\.")
    spdi_re = re.compile(r"(?P<ac>[^:]+):(?P<pos>\d+):(?P<del_len_or_seq>\w+):(?P<ins_seq>\w+)")


    def __init__(self,
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


    def translate_from(self, var, fmt=None):
        """Translate variation `var` to VRS object

        If `fmt` is None, guess the appropriate format and return the variant.
        If `fmt` is specified, try only that format.

        See also notes about `from_` and `to_` methods.
        """

        if fmt:
            t = self.from_translators[fmt]
            o = t(self, var)
            if o is None:
                raise ValueError(f"Unable to parse data as {fmt} variation")
            return o

        for fmt, t in self.from_translators.items():
            o = t(self, var)
            if o:
                return o

        formats = list(self.from_translators.keys())
        raise ValueError(f"Unable to parse data as {', '.join(formats)}")

    def translate_to(self, vo, fmt):
        """translate vrs object `vo` to named format `fmt`"""
        t = self.to_translators[fmt]
        return t(self, vo)


    ############################################################################
    ## INTERNAL

    def _from_beacon(self, beacon_expr, assembly_name=None):
        """Parse beacon expression into VRS Allele

        #>>> a = tlr.from_beacon("13 : 32936732 G > C")
        #>>> a.as_dict()
        {'location': {'interval': {
           'end': {'value': 32936732, 'type': Number},
           'start': {'value': 32936731, 'type': Number},
           'type': 'SequenceInterval'},
          'sequence_id': 'GRCh38:13 ',
          'type': 'SequenceLocation'},
         'state': {'sequence': 'C', 'type': 'LiteralSequenceExpression'},
         'type': 'Allele'}

        """

        if not isinstance(beacon_expr, str):
            return None
        m = self.beacon_re.match(beacon_expr.replace(" ", ""))
        if not m:
            return None

        g = m.groupdict()
        if assembly_name is None:
            assembly_name = self.default_assembly_name
        sequence_id = assembly_name + ":" + g["chr"]
        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        interval = models.SequenceInterval(start=models.Number(value=start),
                                           end=models.Number(value=end))
        location = models.Location(sequence_id=sequence_id, interval=interval)
        sstate = models.LiteralSequenceExpression(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_gnomad(self, gnomad_expr, assembly_name=None):
        """Parse gnomAD-style VCF expression into VRS Allele

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
        if assembly_name is None:
            assembly_name = self.default_assembly_name
        sequence_id = assembly_name + ":" + g["chr"]
        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        interval = models.SequenceInterval(start=models.Number(value=start),
                                           end=models.Number(value=end))
        location = models.Location(sequence_id=sequence_id, interval=interval)
        sstate = models.LiteralSequenceExpression(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_hgvs(self, hgvs_expr):
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

        sv = self._hgvs_parser.parse_hgvs_variant(hgvs_expr)

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

        location = models.Location(sequence_id=sequence_id, interval=interval)
        sstate = models.LiteralSequenceExpression(sequence=state)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_spdi(self, spdi_expr):

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
        location = models.Location(sequence_id=sequence_id, interval=interval)
        sstate = models.LiteralSequenceExpression(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


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

    def _get_hgvs_tools(self):
        """ Only create UTA db connection if needed. There will be one connectionn per translator.
        """
        if self.hgvs_tools is None:
            uta_conn = hgvs.dataproviders.uta.connect()
            self.hgvs_tools = HgvsTools(uta_conn)
        return self.hgvs_tools

    def _to_hgvs(self, vo, namespace="refseq"):
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

        if (type(vo).__name__ != "Allele"
            or type(vo.location).__name__ != "SequenceLocation"
            or type(vo.state).__name__ != "LiteralSequenceExpression"):
            raise ValueError(f"_to_hgvs requires a VRS Allele with SequenceLocation and LiteralSequenceExpression")

        sequence_id = str(vo.location.sequence_id)
        aliases = self.data_proxy.translate_sequence_identifier(sequence_id, namespace)

        # infer type of sequence based on accession
        # TODO: move to bioutils
        stypes = list(set(t for t in (ir_stype(a) for a in aliases) if t))
        if len(stypes) != 1:
            raise ValueError(f"Couldn't infer sequence type for {sequence_id} ({stypes})")
        stype = stypes[0]

        # build interval and edit depending on sequence type
        if stype == "p":
            raise ValueError("Only nucleic acid variation is currently supported")
            # ival = hgvs.location.Interval(start=start, end=end)
            # edit = hgvs.edit.AARefAlt(ref=None, alt=vo.state.sequence)
        else:                   # pylint: disable=no-else-raise
            start = vo.location.interval.start.value
            end = vo.location.interval.end.value
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

            if not (any(a.startswith(pfx) for pfx in ("NM", "NP", "NC", "NG"))):
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
                _logger.warning(f"No data found for accession {a}")

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

        if (type(vo).__name__ != "Allele"
            or type(vo.location).__name__ != "SequenceLocation"
            or type(vo.state).__name__ != "LiteralSequenceExpression"):
            raise ValueError(f"_to_hgvs requires a VRS Allele with SequenceLocation and LiteralSequenceExpression")

        sequence_id = str(vo.location.sequence_id)
        aliases = self.data_proxy.translate_sequence_identifier(sequence_id, namespace)
        aliases = [a.split(":")[1] for a in aliases]

        start = vo.location.interval.start.value
        end = vo.location.interval.end.value
        spdi_tail = f":{start}:{end-start}:{vo.state.sequence}"
        spdis = [a + spdi_tail for a in aliases]
        return spdis


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
            allele = normalize(allele, self.data_proxy)

        if self.identify:
            allele._id = ga4gh_identify(allele)

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
        "13 : 32936732 G > C",
        "NC_000013.11:g.32936732G>C",
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
