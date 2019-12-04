import gzip
import logging
import os

import lxml.etree as le
import smart_open


assert os.environ["USER"] == "reece", "This file is in extras for Reece's convenience. You shouldn't be using it. Go away."


_logger = logging.getLogger()


class VariationArchive:
    def __init__(self, element):
        assert element.tag == "VariationArchive", "Expected node type `VariationArchive`"
        self._element = element

    @property
    def accession(self):
        return self._element.get("Accession")

    @property
    def acv(self):
        return self.accession + "." + self.version

    @property
    def hgvs_expressions(self):
        return [e for e in self._element.xpath(".//Expression/text()")]

    @property
    def is_current(self):
        return self.record_status == "current"

    @property
    def record_status(self):
        return self._element.xpath("RecordStatus/text()")[0]

    @property
    def record_type(self):
        return self._element.get("RecordType")

    @property
    def variation_name(self):
        return self._element.get("VariationName")
    
    @property
    def variation_type(self):
        return self._element.get("VariationType")

    @property
    def version(self):
        return self._element.get("Version")
    
    @property
    def xrefs(self):
        return [e.attrib for e in self._element.xpath(".//XRefList/XRef")]

    def __del__(self):
        self._element.clear()


class ClinvarParser:
    def __init__(self, fp, current_only=True):
        self._xp = le.iterparse(gzip.open(fp), tag="VariationArchive")
        self.current_only = current_only

    def __iter__(self):
        yield from (va for va in (VariationArchive(e) for _, e in self._xp) if not self.current_only or va.is_current)


def clinvar_open(fp):
    yield from ClinvarParser(fp)



if __name__ == "__main__":
    import sys
    fn = sys.argv[1]
    for va in clinvar_open(fn):
        print(f"{va.acv}\t{va.variation_name}\t{va.variation_type}\t{len(va.hgvs_expressions)}\t{len(va.xrefs)}")
