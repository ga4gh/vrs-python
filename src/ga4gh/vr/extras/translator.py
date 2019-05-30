"""Translates various formats into VMC models.

"""

import logging
import re

from bioutils.accessions import coerce_namespace
import hgvs.parser

from ga4gh.vr import models, identify
from .decorators import lazy_property

_logger = logging.getLogger(__name__)


beacon_re = re.compile(r"(?P<chr>[^-]+):(?P<pos>\d+)(?P<ref>\w+)>(?P<alt>\w+)")
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
        
        >>> a = tlr.from_beacon("13 : 32936732 G > C")
        >>> a.as_dict()
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
        allele.id = identify(allele)
        return allele


    def from_hgvs(self, hgvs_expr):
        """parse hgvs into a VR object (typically an Allele)

        >>> a = tlr.from_hgvs("NM_012345.6:c.22A>T")
        >>> a.as_dict()
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
        allele.id = identify(allele)
        return allele

    
    def from_spdi(self, spdi_expr):
        """Parse SPDI expression in to a GA4GH Allele

        >>> a = tlr.from_spdi("NM_012345.6:21:1:T")
        >>> a.as_dict()
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
        allele.id = identify(allele)
        return allele


    def from_vcf(self, vcf_expr, assembly_name=None):
        """Parse gnomAD-style VCF expression into VR Allele

        >>> a = tlr.from_vcf("1-55516888-G-GA")
        >>> a.as_dict()
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
        allele.id = identify(allele)
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
            return self.data_proxy.translate_sequence_identifier(ir)
        return ir


if __name__ == "__main__":
    logging.basicConfig(level="DEBUG")

    local = True
    if local:
        from biocommons.seqrepo import SeqRepo
        from ga4gh.vr.extras.dataproxy import SeqRepoDataProxy
        sr = SeqRepo(root_dir="/usr/local/share/seqrepo/latest/")
        dp = SeqRepoDataProxy(sr=sr)
    else:
        from ga4gh.vr.extras.dataproxy import SeqRepoRESTDataProxy
        dp = SeqRepoRESTDataProxy(base_url="http://localhost:5000/seqrepo")

    tlr = Translator(data_proxy=dp)
