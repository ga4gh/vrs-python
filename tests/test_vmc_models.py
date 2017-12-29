import datetime
import json

from vmc import models, computed_id, serialize

i = models.Interval(start=42, end=42)
assert '{"end": 42, "start": 42}' == str(i.serialize())

l = models.Location(sequence_id="VMC:GS_01234", interval=i)
l.id = computed_id(l)
locations = {l.id: l.as_dict()}
assert '{"id": "VMC:GL__e3BmXQJ0CYG0ZtXeZDBzBQkvnetKxgb", "interval": {"end": 42, "start": 42}, "sequence_id": "VMC:GS_01234"}' == str(
    l.serialize())

a = models.Allele(location_id=l.id, state="A")
a.id = computed_id(a)
alleles = {a.id: a.as_dict()}
assert '{"id": "VMC:GA_6DSejQb2b3Of9N5QXSULOTfU_7i9RA-F", "location_id": "VMC:GL__e3BmXQJ0CYG0ZtXeZDBzBQkvnetKxgb", "state": "A"}' == str(
    a.serialize())

h1 = models.Haplotype(allele_ids=[a.id], completeness="PARTIAL")
h1.id = computed_id(h1)
h2 = models.Haplotype(allele_ids=[a.id], completeness="COMPLETE")
h2.id = computed_id(h2)
haplotypes = {h.id: h.as_dict() for h in [h1, h2]}
assert '{"allele_ids": ["VMC:GA_6DSejQb2b3Of9N5QXSULOTfU_7i9RA-F"], "completeness": "PARTIAL", "id": "VMC:GH_dmPKNQRhrnF-XRsqder9RTPz8ZjRtVLg"}' == str(
    h1.serialize())
assert '{"allele_ids": ["VMC:GA_6DSejQb2b3Of9N5QXSULOTfU_7i9RA-F"], "completeness": "COMPLETE", "id": "VMC:GH_p1wCN7-eVQklcv6tg5d2DGVHCX3eCP2p"}' == str(
    h2.serialize())

g = models.Genotype(haplotype_ids=[h1.id, h2.id], completeness="COMPLETE")
g.id = computed_id(g)
genotypes = {g.id: g.as_dict()}
assert '{"completeness": "COMPLETE", "haplotype_ids": ["VMC:GH_dmPKNQRhrnF-XRsqder9RTPz8ZjRtVLg", "VMC:GH_p1wCN7-eVQklcv6tg5d2DGVHCX3eCP2p"], "id": "VMC:GG_9uOweOFteefYsz-Vno_Yaw1itsjg8MI7"}' == g.serialize(
)

identifiers = {
    "VMC:GS_01234": [
        models.Identifier(namespace="NCBI", accession="NM_0123.4"),
        models.Identifier(namespace="VMC", accession="GS_01234")
    ],
    "VMC:GL__e3BmXQJ0CYG0ZtXeZDBzBQkvnetKxgb_": [
        models.Identifier(namespace="dbSNP", accession="rs12345"),
    ]
}

bundle = models.Vmcbundle(
    meta=models.Meta(
        generated_at=datetime.datetime.isoformat(datetime.datetime.now()),
        vmc_version=0,
    ),
    locations=locations,
    alleles=alleles,
    haplotypes=haplotypes,
    genotypes=genotypes,
    identifiers=identifiers)

bundle_fn = __file__.replace(".py", ".json")
saved_bundle = models.Vmcbundle(**json.load(open(bundle_fn)))

# fudge the timestamp so that they compare equal
saved_bundle.meta.generated_at = bundle.meta.generated_at

assert bundle == saved_bundle
