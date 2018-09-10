from vmc import models
from vmc.conversions import from_hgvs


class BundleManager(models.Vmcbundle):
    def __init__(self):
        super(models.Vmcbundle, self).__init__(
            alleles = {},
            haplotypes = {},
            genotypes = {},
            identifiers = {},
            meta = {},
            )

    def merge_bundle(self, other):
        self.alleles.update(other.alleles)
        self.haplotypes.update(other.haplotypes)
        self.genotypes.update(other.genotypes)
        self.identifiers.update(other.identifiers)
        self.meta.update(other.meta)

    def add_hgvs(self, hgvs):
        self.merge_bundle(from_hgvs(hgvs))

    def verify(self):
        pass



if __name__ == "__main__":
    bm = BundleManager()
    hb = from_hgvs("NM_000314.4:c.493A>C")
