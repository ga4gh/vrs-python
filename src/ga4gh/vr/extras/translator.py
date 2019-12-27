"""Translates various formats into VR models.

"""

import logging
import re

from bioutils.accessions import coerce_namespace
import hgvs.parser

import hgvs.location
import hgvs.posedit
import hgvs.edit
import hgvs.sequencevariant

from ga4gh.core import ga4gh_identify
from ga4gh.vr import models
from .decorators import lazy_property


_logger = logging.getLogger(__name__)


beacon_re = re.compile(r"(?P<chr>[^-]+)\s*:\s*(?P<pos>\d+)\s*(?P<ref>\w+)\s*>\s*(?P<alt>\w+)")
vcf_re = re.compile(r"(?P<chr>[^-]+)-(?P<pos>\d+)-(?P<ref>\w+)-(?P<alt>\w+)")
spdi_re = re.compile(r"(?P<ac>[^:]+):(?P<pos>\d+):(?P<del_len>\d+):(?P<ins_seq>\w+)")



class Translator:
    """Translates various variation formats to and from GA4GH VR models

    """

    def __init__(self,
                 data_proxy,
                 default_assembly_name="GRCh38",
                 translate_sequence_identifiers=True,
                 identify=True):
        self.default_assembly_name = default_assembly_name
        self.translate_sequence_identifiers = translate_sequence_identifiers
        self.data_proxy = data_proxy
        self.identify = identify


    def from_beacon(self, beacon_expr, assembly_name=None):
        """Parse beacon expression into VR Allele
        
        #>>> a = tlr.from_beacon("13 : 32936732 G > C")
        #>>> a.as_dict()
        {'location': {'interval': {'end': 32936732,
           'start': 32936731,
           'type': 'SimpleInterval'},
          'sequence_id': 'GRCh38:13 ',
          'type': 'SequenceLocation'},
         'state': {'sequence': 'C', 'type': 'SequenceState'},
         'type': 'Allele'}

        """          
        
        m = beacon_re.match(beacon_expr.replace(" ", ""))
        if not m:
            raise ValueError(f"Not a Beacon expression: {beacon_expr}")

        g = m.groupdict()
        if assembly_name is None:
            assembly_name = self.default_assembly_name
        sequence_id = assembly_name + ":" + g["chr"]
        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        interval = models.SimpleInterval(start=start, end=end)
        location = models.Location(sequence_id=self._seq_id_mapper(sequence_id), interval=interval)
        sstate = models.SequenceState(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        if self.identify:
            allele._id = ga4gh_identify(allele)
        return allele


    def from_hgvs(self, hgvs_expr):
        """parse hgvs into a VR object (typically an Allele)

        #>>> a = tlr.from_hgvs("NM_012345.6:c.22A>T")
        #>>> a.as_dict()
        {
          'location': {
            'interval': {'end': 22, 'start': 21, 'type': 'SimpleInterval'},
            'sequence_id': 'refseq:NM_012345.6',
            'type': 'SequenceLocation'
          },
          'state': {'sequence': 'T', 'type': 'SequenceState'},
          'type': 'Allele'
        }

        """
        sv = self._hgvs_parser.parse_hgvs_variant(hgvs_expr)

        # prefix accession with namespace
        sequence_id = coerce_namespace(sv.ac)

        if isinstance(sv.posedit.pos, hgvs.location.BaseOffsetInterval):
            if sv.posedit.pos.start.is_intronic or sv.posedit.pos.end.is_intronic:
                raise ValueError("Intronic HGVS variants are not supported ({sv.posedit})")

        if sv.posedit.edit.type == 'ins':
            interval = models.SimpleInterval(start=sv.posedit.pos.start.base,
                                        end=sv.posedit.pos.start.base)
            state = sv.posedit.edit.alt
        elif sv.posedit.edit.type in ('sub', 'del', 'delins', 'identity'):
            interval = models.SimpleInterval(start=sv.posedit.pos.start.base - 1,
                                       end=sv.posedit.pos.end.base)
            if sv.posedit.edit.type == 'identity':
                state = self.data_proxy.get_sequence(sv.ac,
                                                     sv.posedit.pos.start.base - 1,
                                                     sv.posedit.pos.end.base)
            else:
                state = sv.posedit.edit.alt or ''
        else:
            raise ValueError(f"HGVS variant type {sv.posedit.edit.type} is unsupported")

        location = models.Location(sequence_id=self._seq_id_mapper(sequence_id), interval=interval)
        sstate = models.SequenceState(sequence=state)
        allele = models.Allele(location=location, state=sstate)
        if self.identify:
            allele._id = ga4gh_identify(allele)
        return allele

    
    def to_hgvs(self, allele, namespace=None):
        """generates a list of HGVS expressions for VR Allele.

        If `namespace` is not None, returns HGVS strings for the
        specified namespace.

        If `namespace` is None, returns HGVS strings for all alias
        translations.

        If no alias translations are available, an empty list is
        returned.

        If the VR object cannot be expressed as HGVS, raises ValueError.

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


        if (type(allele).__name__ != "Allele"
            or type(allele.location).__name__ != "SequenceLocation"
            or type(allele.state).__name__ != "SequenceState"):
            raise ValueError(f"to_hgvs requires a VR Allele with SequenceLocation and SequenceState")

        sequence_id = str(allele.location.sequence_id)
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
            # edit = hgvs.edit.AARefAlt(ref=None, alt=allele.state.sequence)
        else:
            start = allele.location.interval.start
            end = allele.location.interval.end
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
            alt = str(allele.state.sequence) or None  # "" => None
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
            if ns.startswith("GRC") and namespace is None:
                continue
            var.ac = a
            hgvs_exprs += [str(var)]

        return list(set(hgvs_exprs))

    
    def from_spdi(self, spdi_expr):

        """Parse SPDI expression in to a GA4GH Allele

        #>>> a = tlr.from_spdi("NM_012345.6:21:1:T")
        #>>> a.as_dict()
        {
          'location': {
            'interval': {'end': 22, 'start': 21, 'type': 'SimpleInterval'},
            'sequence_id': 'refseq:NM_012345.6',
            'type': 'SequenceLocation'
          },
          'state': {'sequence': 'T', 'type': 'SequenceState'},
          'type': 'Allele'
        }
        """

        m = spdi_re.match(spdi_expr)
        if not m:
            raise ValueError(f"Not a SPDI expression: {spdi_expr}")

        g = m.groupdict()
        sequence_id = coerce_namespace(g["ac"])
        start = int(g["pos"])
        del_len = int(g["del_len"])
        end = start + del_len
        ins_seq = g["ins_seq"]

        interval = models.SimpleInterval(start=start, end=end)
        location = models.Location(sequence_id=self._seq_id_mapper(sequence_id), interval=interval)
        sstate = models.SequenceState(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        if self.identify:
            allele._id = ga4gh_identify(allele)
        return allele


    def from_vcf(self, vcf_expr, assembly_name=None):
        """Parse gnomAD-style VCF expression into VR Allele

        #>>> a = tlr.from_vcf("1-55516888-G-GA")
        #>>> a.as_dict()
        {'location': {'interval': {'end': 55516888,
           'start': 55516887,
           'type': 'SimpleInterval'},
          'sequence_id': 'GRCh38:1',
          'type': 'SequenceLocation'},
         'state': {'sequence': 'GA', 'type': 'SequenceState'},
         'type': 'Allele'}
        
        """          
        
        m = vcf_re.match(vcf_expr)
        if not m:
            raise ValueError(f"Not a vcf/gnomAD expression: {vcf_expr}")

        g = m.groupdict()
        if assembly_name is None:
            assembly_name = self.default_assembly_name
        sequence_id = assembly_name + ":" + g["chr"]
        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        interval = models.SimpleInterval(start=start, end=end)
        location = models.Location(sequence_id=self._seq_id_mapper(sequence_id), interval=interval)
        sstate = models.SequenceState(sequence=ins_seq)
        allele = models.Allele(location=location, state=sstate)
        if self.identify:
            allele._id = ga4gh_identify(allele)
        return allele


    
    ############################################################################
    ## INTERNAL

    @lazy_property
    def _hgvs_parser(self):
        """instantiates and returns an hgvs parser instance"""
        _logger.info("Creating  parser")
        return hgvs.parser.Parser()


    def _seq_id_mapper(self, ir):
        if self.translate_sequence_identifiers:
            return self.data_proxy.translate_sequence_identifier(ir, "ga4gh")[0]
        return ir
