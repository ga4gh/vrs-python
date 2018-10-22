import datetime
import json

from vmc import models, computed_id, serialize


# Interval
i = models.Interval(start=42, end=42)
assert "<Interval|42|42>" == serialize(i)
assert {"end": 42, "start": 42} == i.as_dict()


# Location
l = models.Location(sequence_id="VMC:GS_01234", interval=i)
assert "<Location|VMC:GS_01234|<Interval|42|42>>" == serialize(l)
l.id = computed_id(l)
assert "VMC:GL_OUqODzxryILUEDmv7uF8R8NwREJAx7gN" == l.id
assert {"id": "VMC:GL_OUqODzxryILUEDmv7uF8R8NwREJAx7gN", "interval": {"end": 42, "start": 42}, "sequence_id": "VMC:GS_01234"} == l.as_dict()

locations = {l.id: l.as_dict()}

# Allele
a = models.Allele(location_id=l.id, state="A")
assert "<Allele|VMC:GL_OUqODzxryILUEDmv7uF8R8NwREJAx7gN|A>" == serialize(a)
a.id = computed_id(a)
assert "VMC:GA_xTR0mmMviMLoAI9SwmDMFYr_AZczkjyU" == a.id
assert {'id': 'VMC:GA_xTR0mmMviMLoAI9SwmDMFYr_AZczkjyU', 'location_id': 'VMC:GL_OUqODzxryILUEDmv7uF8R8NwREJAx7gN', 'state': 'A'} == a.as_dict()

alleles = {a.id: a.as_dict()}


# Haplotype
h1 = models.Haplotype(allele_ids=[a.id], completeness="PARTIAL")
assert "<Haplotype||PARTIAL|[VMC:GA_xTR0mmMviMLoAI9SwmDMFYr_AZczkjyU]>" == serialize(h1)
h1.id = computed_id(h1)
assert "VMC:GH_RhBraQWtIyTu9VP90IJ2wvqOk4tN76SB" == h1.id
assert {"allele_ids": ["VMC:GA_xTR0mmMviMLoAI9SwmDMFYr_AZczkjyU"], "completeness": "PARTIAL", "id": "VMC:GH_RhBraQWtIyTu9VP90IJ2wvqOk4tN76SB"} == h1.as_dict()

h2 = models.Haplotype(allele_ids=[a.id], completeness="COMPLETE")
assert "<Haplotype||COMPLETE|[VMC:GA_xTR0mmMviMLoAI9SwmDMFYr_AZczkjyU]>" == serialize(h2)
h2.id = computed_id(h2)
assert "VMC:GH_x29VkAQwgr724e5kjgylAoWQfk2LqiHd" == h2.id
assert {"allele_ids": ["VMC:GA_xTR0mmMviMLoAI9SwmDMFYr_AZczkjyU"], "completeness": "COMPLETE", "id": "VMC:GH_x29VkAQwgr724e5kjgylAoWQfk2LqiHd"} == h2.as_dict()

haplotypes = {h.id: h.as_dict() for h in [h1, h2]}


# Genotype
g = models.Genotype(haplotype_ids=[h1.id, h2.id], completeness="COMPLETE")
assert "<Genotype|COMPLETE|[VMC:GH_RhBraQWtIyTu9VP90IJ2wvqOk4tN76SB;VMC:GH_x29VkAQwgr724e5kjgylAoWQfk2LqiHd]>" == serialize(g)
g.id = computed_id(g)
assert "VMC:GG_fPYzC18fCsTBXfG48apLSvr2TeEepTRB" == g.id
assert {"completeness": "COMPLETE", "haplotype_ids": ["VMC:GH_RhBraQWtIyTu9VP90IJ2wvqOk4tN76SB", "VMC:GH_x29VkAQwgr724e5kjgylAoWQfk2LqiHd"], "id": "VMC:GG_fPYzC18fCsTBXfG48apLSvr2TeEepTRB"} == g.as_dict()

genotypes = {g.id: g.as_dict()}


# Identifier
identifiers = {
    "VMC:GS_01234": ["RefSeq:NM_0123.4", "VMC:GS_01234"],
    "VMC:GL__e3BmXQJ0CYG0ZtXeZDBzBQkvnetKxgb_": ["dbSNP:rs12345", "VMC:GL__e3BmXQJ0CYG0ZtXeZDBzBQkvnetKxgb_"],
    }


# Bundle
bundle = models.Vmcbundle(
    meta=models.Meta(
        generated_at=datetime.datetime.isoformat(datetime.datetime.now()),
        vmc_version=0,
    ),
    locations=locations,
    alleles=alleles,
    haplotypes=haplotypes,
    genotypes=genotypes,
    identifiers=identifiers,
    )

bundle_fn = __file__.replace(".py", ".json")
# to create (py 3 req'd):
# json.dump(bundle.for_json(), open(bundle_fn + "-new", "w"), indent=3, sort_keys=True)
saved_bundle = models.Vmcbundle(**json.load(open(bundle_fn)))

# fudge the timestamp so that they compare equal
saved_bundle.meta.generated_at = bundle.meta.generated_at

bundle_d = bundle.as_dict()
saved_bundle_d = saved_bundle.as_dict()
assert bundle_d == saved_bundle_d
