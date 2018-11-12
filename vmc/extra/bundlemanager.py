"""manages a VMC bundle and provides helpful a helpful interface to
it"""

import collections
import datetime
import functools
import uuid

import hgvs
import hgvs.edit
import hgvs.posedit
import hgvs.parser
import hgvs.location
import hgvs.sequencevariant

from vmc import models, computed_id, get_vmc_sequence_identifier
from vmc.extra.seqrepo import get_reference_sequence


# TODO: Implement changeable id style: vmcdigest, serial, uuid


# TODO: Use new `from hgvs.easy import parser`
hp = None


def _get_hgvs_parser():
    global hp
    if hp is None:
        hp = hgvs.parser.Parser()
    return hp


def ns_to_name(ns):
    # totally hokey
    if ns == "refseq":
        return "RefSeq"
    return ns



_object_id = 0
def _get_id_serial(o):
    global _object_id
    _object_id += 1
    return str(_object_id)

_id_functions = {
    'computed': computed_id,
    'serial': _get_id_serial,
    'uuid': lambda: str(uuid.uuid4()),
    }



class BundleManager:
    def __init__(self, id_function="computed"):
        self.alleles = {}
        self.genotypes = {}
        self.haplotypes = {}
        self.identifiers = collections.defaultdict(set)
        self.locations = {}
        self._id_function = _id_functions[id_function]


    def add_hgvs_allele(self, hgvs_allele):
        """parse and add the hgvs_allele to the bundle"""
        hp = _get_hgvs_parser()

        sv = hp.parse_hgvs_variant(hgvs_allele)

        sequence_id = get_vmc_sequence_identifier(sv.ac)
        self.identifiers[sequence_id].add(sv.ac)

        if isinstance(sv.posedit.pos, hgvs.location.BaseOffsetInterval):
            if sv.posedit.pos.start.is_intronic or sv.posedit.pos.end.is_intronic:
                raise ValueError("Intronic HGVS variants are not supported".format(sv.posedit.edit.type))

        if sv.posedit.edit.type == 'ins':
            interval = models.Interval(start=sv.posedit.pos.start.base, end=sv.posedit.pos.start.base)
            state = sv.posedit.edit.alt
        elif sv.posedit.edit.type in ('sub', 'del', 'delins', 'identity'):
            interval = models.Interval(start=sv.posedit.pos.start.base - 1, end=sv.posedit.pos.end.base)
            if sv.posedit.edit.type == 'identity':
                state = get_reference_sequence(sv.ac, sv.posedit.pos.start.base - 1, sv.posedit.pos.end.base)
            else:
                state = sv.posedit.edit.alt or ''
        else:
            raise ValueError("HGVS variant type {} is unsupported".format(sv.posedit.edit.type))

        location = models.Location(sequence_id=sequence_id, interval=interval)
        location.id = self._id_function(location)
        self.locations[location.id] = location

        allele = models.Allele(location_id=location.id, state=state)
        allele.id = self._id_function(allele)
        self.alleles[allele.id] = allele

        return allele


    def add_hgvs_haplotype(self, hgvs_alleles, completeness="UNKNOWN"):
        alleles = [self.add_hgvs_allele(hgvs_allele) for hgvs_allele in hgvs_alleles]

        # create location from bounding box around alleles
        sequence_ids = set(self.locations[a.location_id].sequence_id for a in alleles)
        if len(sequence_ids) > 1:
            raise Exception("Haplotypes must be defined on a single sequence")
        sequence_id = next(iter(sequence_ids))
        intervals = [self.locations[a.location_id].interval for a in alleles]
        interval_min = min(int(i.start) for i in intervals)
        interval_max = max(int(i.end) for i in intervals)
        interval = models.Interval(start=interval_min, end=interval_max)
        location = models.Location(sequence_id=sequence_id, interval=interval)
        location.id = self._id_function(location)
        self.locations[location.id] = location

        haplotype = models.Haplotype(completeness=completeness,
                                     location_id=location.id,
                                     allele_ids=[a.id for a in alleles])
        haplotype.id = self._id_function(haplotype)
        self.haplotypes[haplotype.id] = haplotype
        return haplotype


    def add_hgvs_genotype(self, hgvs_haplotypes, completeness="UNKNOWN"):
        haplotypes = [self.add_hgvs_haplotype(hh) for hh in hgvs_haplotypes]
        genotype = models.Genotype(completeness=completeness,
                                   haplotype_ids=[h.id for h in haplotypes])
        genotype.id = self._id_function(genotype)
        self.genotypes[genotype.id] = genotype
        return genotype


    def as_hgvs(self):
        """returns a list of HGVS alleles, haplotypes, and genotypes"""

        @functools.lru_cache()
        def allele_as_hgvs(allele_id):
            allele = self.alleles[allele_id]
            location = self.locations[allele.location_id]
            interval = location.interval
            seq_ir = next(iter(self.identifiers[location.sequence_id]))
            ns, acc = seq_ir.split(":") if ":" in seq_ir else (None, seq_ir)
            ref = get_reference_sequence(seq_ir, interval.start, interval.end)
            type = "g"
            v = hgvs.sequencevariant.SequenceVariant(
                ac = acc,
                type = type,
                posedit = hgvs.posedit.PosEdit(
                    pos = hgvs.location.Interval(
                        start=hgvs.location.SimplePosition(interval.start+1),
                        end=hgvs.location.SimplePosition(interval.end)),
                    edit = hgvs.edit.NARefAlt(ref=ref, alt=allele.state)
                    )
                )
            return str(v)

        return {"alleles": [allele_as_hgvs(aid)
                            for aid in self.alleles.keys()],
                "haplotypes": [[allele_as_hgvs(aid)
                                for aid in h.allele_ids]
                               for h in self.haplotypes.values()],
                "genotypes": [[[allele_as_hgvs(aid)
                                for aid in self.haplotypes[hid].allele_ids]
                               for hid in g.haplotype_ids]
                              for g in self.genotypes.values()]
        }


    def as_bundle(self):
        b = models.Vmcbundle(
            alleles = self.alleles,
            genotypes = self.genotypes,
            haplotypes = self.haplotypes,
            identifiers = {k: list(v) for k, v in self.identifiers.items()},
            locations = self.locations,
            meta=models.Meta(
                generated_at=datetime.datetime.isoformat(datetime.datetime.now()),
                vmc_version=0,
            ))
        return b

    def as_dict(self):
        return self.as_bundle().as_dict()

    def get_referenced_sequence_ids(self):
        return set(l.sequence_id for l in self.locations.values())


    def merge_bundle(self, other):
        self.alleles.update(other.alleles)
        self.genotypes.update(other.genotypes)
        self.haplotypes.update(other.haplotypes)
        self.identifiers.update(other.identifiers)
        self.locations.update(other.locations)
        self.meta.update(other.meta)



if __name__ == "__main__":
    bm = BundleManager()

    allele = bm.add_hgvs_allele("NM_000041.3:c.388T>C")

    haplotype = bm.add_hgvs_haplotype(["NM_000041.3:c.388T>C",
                                       "NM_000041.3:c.526C>T"])
    
    genotype = bm.add_hgvs_genotype([["NM_000041.3:c.388T>C",
                                      "NM_000041.3:c.526C>T"],
                                     ["NM_000041.3:c.388T>C",
                                      "NM_000041.3:c.526C>T"]])

