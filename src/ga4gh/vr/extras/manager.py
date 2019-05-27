"""manages a GA4GH VR bundle and provides a helpful interface to it.

"""

import logging

from ga4gh.vr import computed_id


_logger = logging.getLogger(__name__)


class Manager:
    def __init__(self, storage):
        self.storage = storage

    def add_allele(self, allele):
        allele.id = ga4gh.vr.computed_id(allele)
        self.storage.alleles[allele.id] = allele
        self.add_location(allele.location)
        _logger.warn(f"Added Allele {allele.id}")
        return allele

    def add_location(self, location):
        location.id = ga4gh.vr.computed_id(location)
        self.storage.locations[location.id] = location
        return location

    def add_text(self, text):
        text.id = ga4gh.vr.computed_id(text)
        self.storage.texts[text.id] = text
        return text


    # TODO: replace getters with direct access
    def get_allele(self, id):
        return self.storage.alleles[id]

    def get_location(self, id):
        return self.storage.location[id]

    def get_text(self, id):
        return self.storage.texts[id]
