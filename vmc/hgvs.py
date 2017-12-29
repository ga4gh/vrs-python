from . import models, computed_id
from .seqrepo import get_vmc_sequence_id

import hgvs
import hgvs.parser
import hgvs.location

hp = None


def _get_hgvs_parser():
    global hp
    if hp is None:
        hp = hgvs.parser.Parser()
    return hp


def from_hgvs(hgvs_string):
    hp = _get_hgvs_parser()
    sv = hp.parse_hgvs_variant(hgvs_string)

    ir = models.Identifier(namespace="NCBI", accession=sv.ac)
    sequence_id = get_vmc_sequence_id(ir)

    if isinstance(sv.posedit.pos, hgvs.location.BaseOffsetInterval):
        if sv.posedit.pos.start.is_intronic() or sv.posedit.pos.end.is_intronic():
            raise ValueError("Intronic HGVS variants are not supported".format(sv.posedit.edit.type))
            
    if sv.posedit.edit.type == 'ins':
        interval = models.Interval(start=sv.posedit.pos.start.base,
                                   end=sv.posedit.pos.start.base)
    elif sv.posedit.edit.type in ('sub', 'del', 'delins'):
        interval = models.Interval(start=sv.posedit.pos.start.base - 1,
                                   end=sv.posedit.pos.end.base)
    else:
        raise ValueError("HGVS variant type {} is unsupported".format(sv.posedit.edit.type))
    
    location = models.Location(
        sequence_id = sequence_id,
        interval = interval)
    location.id = computed_id(location)

    allele = models.Allele(
        location_id = location.id,
        state = sv.posedit.edit.alt
        )
    allele.id = computed_id(allele)

    bundle = models.Vmcbundle(
        locations = {location.id: location.as_dict()},
        alleles = {allele.id: allele.as_dict()},
        identifiers = {sequence_id: [ir.as_dict()]},
        )

    return bundle
