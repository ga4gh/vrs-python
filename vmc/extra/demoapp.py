"""Demostrates storing application data adjacent to VMC bundles"""

from collections import defaultdict

from .bundlemanager import BundleManager


class DemoApp:
    def __init__(self):
        self.bm = BundleManager()
        self.annotations = defaultdict(lambda: dict())

    def as_dict(self):
        return {
            "appinfo": {"what": "ever"},
            "vmcbundle": self.bm.as_dict(),
            "annotations": self.annotations,
        }
    
    def get_record(self, id):
        """given variation id, return dictionary of all associated data"""
        def get_object(id):
            if id in self.bm.alleles:
                return self.bm.alleles[id]
            if id in self.bm.haplotypes:
                return self.bm.haplotypes[id]
            if id in self.bm.genotypes:
                return self.bm.haplotypes[id]
            if id in self.bm.locations:
                return self.bm.locations[id]
            return None

        def infer_type(o):
            return type(o).__name__

        o = get_object(id)
        rv = {
            "object": o,
            "type": infer_type(o),
            }
        annotations = {anntype: annots[id] for anntype, annots in self.annotations.items()}
        rv.update(annotations)
        return rv
