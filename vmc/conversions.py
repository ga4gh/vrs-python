# -*- coding: utf-8 -*-

"""generate VMC Bundle from HGVS string

>>> vb = from_hgvs("NC_000019.10:g.44908684C>T")
>>> list(vb.locations.keys())[0]
'VMC:GL_L1IS6jOwSUsOpKihGRcqxHul1IwbV-1s'
>>> list(vb.alleles.keys())[0]
'VMC:GA_AnJl99FJB5tNPupduz8I4R8CCuwCpIY0'

# These types work:
>>> _ = from_hgvs("NM_000314.4:c.706G>T")
>>> _ = from_hgvs("NM_000314.4:c.706_706delG")
>>> _ = from_hgvs("NM_000314.4:c.706_707insTT")
>>> _ = from_hgvs("NM_000314.4:c.706_708delGAC")
>>> _ = from_hgvs("NM_000314.4:c.706_708delGACinsTTGT")
>>> _ = from_hgvs("NM_000314.4:c.706_708delinsTTGT")
>>> _ = from_hgvs("NM_000314.4:c.706delG")

# Some HGVS expressions are not supported:
>>> from_hgvs("NM_000314.4:c.493-2A>C")
Traceback (most recent call last):
...
ValueError: Intronic HGVS variants are not supported

>>> from_hgvs("NM_000314.4:c.493dup")
Traceback (most recent call last):
...
ValueError: HGVS variant type dup is unsupported

"""

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
        meta={"version": "0.1"},
    )

    return bundle
