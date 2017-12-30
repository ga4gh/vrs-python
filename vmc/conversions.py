# -*- coding: utf-8 -*-

"""generate VMC Bundle from HGVS string

>>> vb = from_hgvs("NC_000019.10:g.44908684C>T")
>>> list(vb.locations.keys())[0]
'VMC:GL_9Jht-lguk_jnBvG-wLJbjmBw5v_v7rQo'
>>> list(vb.alleles.values())[0]
<Allele id=VMC:GA_xXBYkzzu1AH0HRbLeFESvllmAKUNN1MF location_id=VMC:GL_9Jht-lguk_jnBvG-wLJbjmBw5v_v7rQo state=T>

>>> vb = from_hgvs("NM_000314.4:c.706_707insTT")

>>> from_hgvs("NM_000314.4:c.493-2A>C")
Traceback (most recent call last):
...
ValueError: Intronic HGVS variants are not supported

>>> from_hgvs("NM_000314.4:c.493dup")
Traceback (most recent call last):
...
ValueError: HGVS variant type dup is unsupported

"""

# vb = from_hgvs("NM_000314.4:c.706G>T")
# vb = from_hgvs("NM_000314.4:c.706_706delG")
# vb = from_hgvs("NM_000314.4:c.706_708delGAC")
# vb = from_hgvs("NM_000314.4:c.706_708delGACinsTTGT")
# vb = from_hgvs("NM_000314.4:c.706_708delinsTTGT")
# vb = from_hgvs("NM_000314.4:c.706delG")
# vb = from_hgvs("NM_000314.4:c.493-2A>C")





import hgvs
import hgvs.parser
import hgvs.location

from . import models, computed_id
from .seqrepo import get_vmc_sequence_id


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
        if sv.posedit.pos.start.is_intronic or sv.posedit.pos.end.is_intronic:
            raise ValueError("Intronic HGVS variants are not supported".format(sv.posedit.edit.type))

    if sv.posedit.edit.type == 'ins':
        interval = models.Interval(start=sv.posedit.pos.start.base, end=sv.posedit.pos.start.base)
    elif sv.posedit.edit.type in ('sub', 'del', 'delins'):
        interval = models.Interval(start=sv.posedit.pos.start.base - 1, end=sv.posedit.pos.end.base)
    else:
        raise ValueError("HGVS variant type {} is unsupported".format(sv.posedit.edit.type))

    location = models.Location(sequence_id=sequence_id, interval=interval)
    location.id = computed_id(location)

    state = sv.posedit.edit.alt or ''
    allele = models.Allele(location_id=location.id, state=state)
    allele.id = computed_id(allele)

    bundle = models.Vmcbundle(
        alleles={allele.id: allele.as_dict()},
        genotypes={},
        haplotypes={},
        identifiers={sequence_id: [ir.as_dict()]},
        locations={location.id: location.as_dict()},
        meta={}
    )

    return bundle
