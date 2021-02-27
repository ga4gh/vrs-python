"""debugging handling when an attribute is sometimes referrable

Background:
Previously, some VRS objects had "referrable" attributes, which means
that they could be either the inline object or a reference.  For
example, Allele.location could be either an inline Location or a
CURIE. Things got more complicated when an attribute used oneOf,
because that required checking whether any of oneOf options was
identifiable; if so, we assumed that they were all identifiable. 



"""


from ga4gh.vrs import models, vrs_enref

sl = models.SequenceLocation(
    sequence_id="ga4gh:SQ.abc123",
    interval=models.SimpleInterval(start=20, end=20))
a = models.Allele(
    location=sl,
    state=models.RepeatedSequence(
        sequence=models.LiteralSequence(sequence="CAG"),
        count={"min": 0, "max": 30}))
g = models.Gene(gene_id="ncbigene:1234")
r = models.IntegerRange(min=0, max=5)

for vo in (
        models.AbsoluteAbundance(subject=a, amount=r),
        # models.AbsoluteAbundance(subject=g, amount=r),
        models.AbsoluteAbundance(subject="ncbigene:1234", amount=r),
        ):
    print(vo.as_dict())
    print(vrs_enref(vo).as_dict())
