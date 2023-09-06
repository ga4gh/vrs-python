"""Translates various external formats into VRS models.

Input formats: VRS (serialized), hgvs, spdi, gnomad (vcf), beacon
Output formats: VRS (serialized), hgvs, spdi, gnomad (vcf)

"""

from collections.abc import Mapping
from typing import Union
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
    spdi_re = re.compile(r"(?P<ac>[^:]+):(?P<pos>\d+):(?P<del_len_or_seq>\w*):(?P<ins_seq>\w*)")


    def __init__(self,
                 data_proxy,
                 default_assembly_name="GRCh38",
                 normalize=True,
                 identify=True):
        self.default_assembly_name = default_assembly_name
        self.data_proxy = data_proxy
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
    # INTERNAL

    def _get_refget_accession(self, alias):
        """Get refget accession given alias

        :param alias: Alias to get accession for
        :return: Refget Accession if found
        """
        refget_accession = None
        try:
            aliases = self.data_proxy.translate_sequence_identifier(
                alias, namespace="ga4gh"
            )
        except KeyError:
            pass
        else:
            if aliases:
                refget_accession = aliases[0].split("ga4gh:")[-1]

        return refget_accession

    def _from_beacon(self, beacon_expr, assembly_name=None):
        """Parse beacon expression into VRS Allele

        #>>> a = tlr.from_beacon("13 : 32936732 G > C")
        #>>> a.model_dump()
        {
          'location': {
            'end': 32936732,
            'start': 32936731,,
            'sequenceReference': {
              'type': 'SequenceReference',
              'refgetAccession': 'SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT'
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
        if assembly_name is None:
            assembly_name = self.default_assembly_name
        sequence = assembly_name + ":" + g["chr"]
        refget_accession = self._get_refget_accession(sequence)
        if not refget_accession:
            return None

        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        seq_ref = models.SequenceReference(refgetAccession=refget_accession)
        location = models.SequenceLocation(sequenceReference=seq_ref, start=start, end=end)
        state = models.LiteralSequenceExpression(sequence=ins_seq)
        allele = models.Allele(location=location, state=state)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_gnomad(self, gnomad_expr, assembly_name=None):
        """Parse gnomAD-style VCF expression into VRS Allele

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
        if assembly_name is None:
            assembly_name = self.default_assembly_name
        sequence = assembly_name + ":" + g["chr"]
        refget_accession = self._get_refget_accession(sequence)
        if not refget_accession:
            return None

        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        seq_ref = models.SequenceReference(refgetAccession=refget_accession)
        location = models.SequenceLocation(sequenceReference=seq_ref, start=start, end=end)
        sstate = models.LiteralSequenceExpression(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_hgvs(self, hgvs_expr):
        """parse hgvs into a VRS object (typically an Allele)

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

        if not isinstance(hgvs_expr, str):
            return None
        if not self.hgvs_re.match(hgvs_expr):
            return None

        sv = self._hgvs_parser.parse_hgvs_variant(hgvs_expr)

        # prefix accession with namespace
        sequence = coerce_namespace(sv.ac)
        refget_accession = self._get_refget_accession(sequence)
        if not refget_accession:
            return None

        if isinstance(sv.posedit.pos, hgvs.location.BaseOffsetInterval):
            if sv.posedit.pos.start.is_intronic or sv.posedit.pos.end.is_intronic:
                raise ValueError("Intronic HGVS variants are not supported ({sv.posedit})")

        if sv.posedit.edit.type == "ins":
            start = sv.posedit.pos.start.base
            end = sv.posedit.pos.start.base
            state = sv.posedit.edit.alt
        elif sv.posedit.edit.type in ("sub", "del", "delins", "identity"):
            start = sv.posedit.pos.start.base - 1
            end = sv.posedit.pos.end.base
            if sv.posedit.edit.type == "identity":
                state = self.data_proxy.get_sequence(sv.ac,
                                                     sv.posedit.pos.start.base - 1,
                                                     sv.posedit.pos.end.base)
            else:
                state = sv.posedit.edit.alt or ""
        elif sv.posedit.edit.type == "dup":
            start = sv.posedit.pos.start.base - 1
            end = sv.posedit.pos.end.base
            ref = self.data_proxy.get_sequence(sv.ac,
                                               sv.posedit.pos.start.base - 1,
                                               sv.posedit.pos.end.base)
            state = ref + ref
        else:
            raise ValueError(f"HGVS variant type {sv.posedit.edit.type} is unsupported")

        seq_ref = models.SequenceReference(refgetAccession=refget_accession)
        location = models.SequenceLocation(sequenceReference=seq_ref, start=start, end=end)
        sstate = models.LiteralSequenceExpression(sequence=state)
        allele = models.Allele(location=location, state=sstate)
        allele = self._post_process_imported_allele(allele)
        return allele


    def _from_spdi(self, spdi_expr):

        """Parse SPDI expression in to a GA4GH Allele

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
        sequence = coerce_namespace(g["ac"])
        refget_accession = self._get_refget_accession(sequence)
        if not refget_accession:
            return None

        start = int(g["pos"])
        try:
            del_len = int(g["del_len_or_seq"])
        except ValueError:
            del_len = len(g["del_len_or_seq"])
        end = start + del_len
        ins_seq = g["ins_seq"]

        seq_ref = models.SequenceReference(refgetAccession=refget_accession)
        location = models.SequenceLocation(sequenceReference=seq_ref, start=start, end=end)
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

        This method assumes that IRIs are dereferenced, providing a `SequenceReference`
        as the `vo.location.sequenceReference`. If a `SequenceReference` is not
        provided, raises TypeError
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

        if not isinstance(vo.location.sequenceReference, models.SequenceReference):
            raise TypeError(
                "`vo.location.sequenceReference` expects a `SequenceReference`"
            )

        sequence = f"ga4gh:{export_sequencelocation_sequence_id(vo.location.sequenceReference)}"
        aliases = self.data_proxy.translate_sequence_identifier(sequence, namespace)

        # infer type of sequence based on accession
        # TODO: move to bioutils
        stypes = list(set(t for t in (ir_stype(a) for a in aliases) if t))
        if len(stypes) != 1:
            raise ValueError(f"Couldn't infer sequence type for {sequence} ({stypes})")
        stype = stypes[0]

        # build interval and edit depending on sequence type
        if stype == "p":
            raise ValueError("Only nucleic acid variation is currently supported")
            # ival = hgvs.location.Interval(start=start, end=end)
            # edit = hgvs.edit.AARefAlt(ref=None, alt=vo.state.sequence)
        else:                   # pylint: disable=no-else-raise
            start, end = vo.location.start, vo.location.end
            # ib: 0 1 2 3 4 5
            #  h:  1 2 3 4 5
            if start == end:    # insert: hgvs uses *exclusive coords*
                ref = None
                end += 1
            else:               # else: hgvs uses *inclusive coords*
                ref = self.data_proxy.get_sequence(sequence, start, end)
                start += 1
            ival = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(start),
                end=hgvs.location.SimplePosition(end)
            )
            alt = str(vo.state.sequence.root) or None  # "" => None
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
        sequence = f"ga4gh:{export_sequencelocation_sequence_id(vo.location.sequenceReference)}"
        aliases = self.data_proxy.translate_sequence_identifier(sequence, namespace)
        aliases = [a.split(":")[1] for a in aliases]
        start, end = vo.location.start, vo.location.end
        spdi_tail = f":{start}:{end-start}:{vo.state.sequence.root}"
        spdis = [a + spdi_tail for a in aliases]
        return spdis

    @lazy_property
    def _hgvs_parser(self):
        """instantiates and returns an hgvs parser instance"""
        _logger.info("Creating  parser")
        return hgvs.parser.Parser()


    def _post_process_imported_allele(self, allele):
        """
        Provide common post-processing for imported Alleles IN-PLACE.
        """
        if self.normalize:
            allele = normalize(allele, self.data_proxy)

        if self.identify:
            allele.id = ga4gh_identify(allele)
            allele.location.id = ga4gh_identify(allele.location)

        return allele


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
        # "gnomad": to_gnomad,
    }


def export_sequencelocation_sequence_id(
    location_sequence_reference: Union[models.IRI, models.SequenceReference]
):
    if isinstance(location_sequence_reference, models.IRI):
        return location_sequence_reference.root
    elif isinstance(location_sequence_reference, models.SequenceReference):
        return location_sequence_reference.refgetAccession


if __name__ == "__main__":
    # pylint: disable=ungrouped-imports

    import coloredlogs
    coloredlogs.install(level="INFO")

    from ga4gh.vrs.dataproxy import create_dataproxy
    # dp = create_dataproxy("seqrepo+file:///usr/local/share/seqrepo/latest")
    dp = create_dataproxy("seqrepo + http://localhost:5555/seqrepo")
    tlr = Translator(data_proxy=dp)

    expressions = [
        "bogus",
        "1-55516888-G-GA",
        "13 : 32936732 G > C",
        "NC_000013.11:g.32936732G>C",
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
