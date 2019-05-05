"""Translates various formats into VMC models.

"""

import logging
import re

import hgvs.parser

import vmc

from ..utils import coerce_namespace
from .decorators import lazy_property


_logger = logging.getLogger(__name__)


beacon_re = re.compile(r"(?P<chr>[^-]+)\s*:\s*(?P<pos>\d+)\s*(?P<ref>\w+)\s*>\s*(?P<alt>\w+)")
gnomad_re = re.compile(r"(?P<chr>[^-]+)-(?P<pos>\d+)-(?P<ref>\w+)-(?P<alt>\w+)")
spdi_re = re.compile(r"(?P<ac>[^:]+):(?P<pos>\d+):(?P<del_len>\d+):(?P<ins_seq>\w+)")



class Translator:
    """Translates various variation formats to and from GA4GH VR models

    """

    def __init__(self, default_assembly_name="GRCh38"):
        # TODO: add identify option
        # TODO: add optional ac translation
        self.default_assembly_name = default_assembly_name

    @lazy_property
    def hgvs_parser(self):
        """Parsing the hgvs grammar is slow (~1s). :-("""
        _logger.warning("Creating parser")
        _logger.info("Creating  parser")
        return hgvs.parser.Parser()


    def from_beacon(self, beacon_expr, assembly_name=None):
        """Parse beacon expression into VR Allele
        
        >>> a = tlr.from_beacon("13 : 32936732 G > C")
        >>> a.as_dict()
        {'location': {'region': {'end': 32936732,
           'start': 32936731,
           'type': 'SimpleInterval'},
          'sequence_id': 'GRCh38:13 ',
          'type': 'SequenceLocation'},
         'state': {'sequence': 'C', 'type': 'SequenceState'},
         'type': 'Allele'}

        """          
        
        m = beacon_re.match(beacon_expr)
        if not m:
            raise ValueError(f"Not a gnomAD expression: {gnomad_expr}")

        g = m.groupdict()
        if assembly_name is None:
            assembly_name = self.default_assembly_name
        sequence_id = assembly_name + ":" + g["chr"]
        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        region = ga4gh.vr.models.SimpleInterval(start=start, end=end)
        location = ga4gh.vr.models.Location(sequence_id=sequence_id, region=region)
        sstate = ga4gh.vr.models.SequenceState(sequence=ins_seq)
        allele = ga4gh.vr.models.Allele(location=location, state=sstate)
        
        return allele


    def from_gnomad(self, gnomad_expr, assembly_name=None):
        """Parse gnomAD expression into VR Allele
        
        >>> a = tlr.from_gnomad("1-55516888-G-GA")
        >>> a.as_dict()
        {'location': {'region': {'end': 55516888,
           'start': 55516887,
           'type': 'SimpleInterval'},
          'sequence_id': 'GRCh38:1',
          'type': 'SequenceLocation'},
         'state': {'sequence': 'GA', 'type': 'SequenceState'},
         'type': 'Allele'}

        """          
        
        m = gnomad_re.match(gnomad_expr)
        if not m:
            raise ValueError(f"Not a gnomAD expression: {gnomad_expr}")

        g = m.groupdict()
        if assembly_name is None:
            assembly_name = self.default_assembly_name
        sequence_id = assembly_name + ":" + g["chr"]
        start = int(g["pos"]) - 1
        ref = g["ref"]
        alt = g["alt"]
        end = start + len(ref)
        ins_seq = alt

        region = ga4gh.vr.models.SimpleInterval(start=start, end=end)
        location = ga4gh.vr.models.Location(sequence_id=sequence_id, region=region)
        sstate = ga4gh.vr.models.SequenceState(sequence=ins_seq)
        allele = ga4gh.vr.models.Allele(location=location, state=sstate)
        
        return allele


    def from_hgvs(self, hgvs_expr):
        """parse hgvs into a VR object (typically an Allele)

        >>> a = tlr.from_hgvs("NM_012345.6:c.22A>T")
        >>> a.as_dict()
        {
          'location': {
            'region': {'end': 22, 'start': 21, 'type': 'SimpleInterval'},
            'sequence_id': 'refseq:NM_012345.6',
            'type': 'SequenceLocation'
          },
          'state': {'sequence': 'T', 'type': 'SequenceState'},
          'type': 'Allele'
        }

        """
        sv = self.hgvs_parser.parse_hgvs_variant(hgvs_expr)

        # prefix accession with namespace
        sequence_id = coerce_namespace(sv.ac)

        if isinstance(sv.posedit.pos, hgvs.location.BaseOffsetInterval):
            if sv.posedit.pos.start.is_intronic or sv.posedit.pos.end.is_intronic:
                raise ValueError("Intronic HGVS variants are not supported ({sv.posedit})")

        if sv.posedit.edit.type == 'ins':
            region = ga4gh.vr.models.SimpleInterval(start=sv.posedit.pos.start.base,
                                        end=sv.posedit.pos.start.base)
            state = sv.posedit.edit.alt
        elif sv.posedit.edit.type in ('sub', 'del', 'delins', 'identity'):
            region = ga4gh.vr.models.SimpleInterval(start=sv.posedit.pos.start.base - 1,
                                       end=sv.posedit.pos.end.base)
            if sv.posedit.edit.type == 'identity':
                state = get_reference_sequence(sv.ac, sv.posedit.pos.start.base - 1,
                                               sv.posedit.pos.end.base)
            else:
                state = sv.posedit.edit.alt or ''
        else:
            raise ValueError(f"HGVS variant type {sv.posedit.edit.type} is unsupported")

        location = ga4gh.vr.models.Location(sequence_id=sequence_id, region=region)
        sstate = ga4gh.vr.models.SequenceState(sequence=state)
        allele = ga4gh.vr.models.Allele(location=location, state=sstate)

        return allele

    
    def from_spdi(self, spdi_expr):
        """Parse SPDI expression in to a GA4GH Allele

        >>> a = tlr.from_spdi("NM_012345.6:21:1:T")
        >>> a.as_dict()
        {
          'location': {
            'region': {'end': 22, 'start': 21, 'type': 'SimpleInterval'},
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

        region = ga4gh.vr.models.SimpleInterval(start=start, end=end)
        location = ga4gh.vr.models.Location(sequence_id=sequence_id, region=region)
        sstate = ga4gh.vr.models.SequenceState(sequence=ins_seq)
        allele = ga4gh.vr.models.Allele(location=location, state=sstate)
        
        return allele
