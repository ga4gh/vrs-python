import datetime
import json

from vmc import models, computed_id, vmc_serialize


# Interval
i = models.Interval(start=42, end=42)
assert "<Interval:42:42>" == vmc_serialize(i)
assert {"end": 42, "start": 42} == i.as_dict()


# Location
l = models.Location(sequence_id="VMC:GS_01234", interval=i)
assert "<Location:<Identifier:VMC:GS_01234>:<Interval:42:42>>" == vmc_serialize(l)
l.id = computed_id(l)
assert "VMC:GL__e3BmXQJ0CYG0ZtXeZDBzBQkvnetKxgb" == l.id
assert {"id": "VMC:GL__e3BmXQJ0CYG0ZtXeZDBzBQkvnetKxgb", "interval": {"end": 42, "start": 42}, "sequence_id": "VMC:GS_01234"} == l.as_dict()

locations = {l.id: l.as_dict()}


# Allele
a = models.Allele(location_id=l.id, state="A")
assert "<Allele:<Identifier:VMC:GL__e3BmXQJ0CYG0ZtXeZDBzBQkvnetKxgb>:A>" == vmc_serialize(a)
a.id = computed_id(a)
assert "VMC:GA_6DSejQb2b3Of9N5QXSULOTfU_7i9RA-F" == a.id
assert {"id": "VMC:GA_6DSejQb2b3Of9N5QXSULOTfU_7i9RA-F", "location_id": "VMC:GL__e3BmXQJ0CYG0ZtXeZDBzBQkvnetKxgb", "state": "A"} == a.as_dict()

alleles = {a.id: a.as_dict()}


# Haplotype
h1 = models.Haplotype(allele_ids=[a.id], completeness="PARTIAL")
assert "<Haplotype::PARTIAL:[<Identifier:VMC:GA_6DSejQb2b3Of9N5QXSULOTfU_7i9RA-F>]>" == vmc_serialize(h1)
h1.id = computed_id(h1)
assert "VMC:GH__oZ5-57j6eI2M12RVJG46a3Q1LS-4IuN" == h1.id
assert {"allele_ids": ["VMC:GA_6DSejQb2b3Of9N5QXSULOTfU_7i9RA-F"], "completeness": "PARTIAL", "id": "VMC:GH__oZ5-57j6eI2M12RVJG46a3Q1LS-4IuN"} == h1.as_dict()

h2 = models.Haplotype(allele_ids=[a.id], completeness="COMPLETE")
assert "<Haplotype::COMPLETE:[<Identifier:VMC:GA_6DSejQb2b3Of9N5QXSULOTfU_7i9RA-F>]>" == vmc_serialize(h2)
h2.id = computed_id(h2)
assert "VMC:GH_zrCwcEZkWF4nDxSk5BQjrEwufUxCXEWq" == h2.id
assert {"allele_ids": ["VMC:GA_6DSejQb2b3Of9N5QXSULOTfU_7i9RA-F"], "completeness": "COMPLETE", "id": "VMC:GH_zrCwcEZkWF4nDxSk5BQjrEwufUxCXEWq"} == h2.as_dict()

haplotypes = {h.id: h.as_dict() for h in [h1, h2]}


# Genotype
g = models.Genotype(haplotype_ids=[h1.id, h2.id], completeness="COMPLETE")
assert "<Genotype:COMPLETE:[<Identifier:VMC:GH__oZ5-57j6eI2M12RVJG46a3Q1LS-4IuN>;<Identifier:VMC:GH_zrCwcEZkWF4nDxSk5BQjrEwufUxCXEWq>]>" == vmc_serialize(g)
g.id = computed_id(g)
assert "VMC:GG_wkq5IDZPN1dTD-WOxdl_8iOB5CAUQ8Va" == g.id
assert {"completeness": "COMPLETE", "haplotype_ids": ["VMC:GH__oZ5-57j6eI2M12RVJG46a3Q1LS-4IuN", "VMC:GH_zrCwcEZkWF4nDxSk5BQjrEwufUxCXEWq"], "id": "VMC:GG_wkq5IDZPN1dTD-WOxdl_8iOB5CAUQ8Va"} == g.as_dict()

genotypes = {g.id: g.as_dict()}


# Identifier
identifiers = {
    "VMC:GS_01234": [
        models.Identifier(namespace="NCBI", accession="NM_0123.4"),
        models.Identifier(namespace="VMC", accession="GS_01234")
    ],
    "VMC:GL__e3BmXQJ0CYG0ZtXeZDBzBQkvnetKxgb_": [
        models.Identifier(namespace="dbSNP", accession="rs12345"),
    ]
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
    identifiers=identifiers)

bundle_fn = __file__.replace(".py", ".json")
# to create (py 3 req'd):
# json.dump(bundle.for_json(), open(bundle_fn + "-new", "w"), indent=3, sort_keys=True)
saved_bundle = models.Vmcbundle(**json.load(open(bundle_fn)))

# fudge the timestamp so that they compare equal
saved_bundle.meta.generated_at = bundle.meta.generated_at

assert bundle == saved_bundle
