import collections
import hashlib
import itertools
import json
import logging
import shelve
import tempfile
import time

import requests

from vmc import models, computed_id, get_vmc_sequence_identifier
from vmc.extra.bundlemanager import BundleManager

_logger = logging.getLogger(__name__)


class RefSeqMapper(shelve.DbfilenameShelf):
    """ClinGen Allele Registry invented their own reference sequence ids,
    which appear to just be synonyms for more established names.  This
    class provides a persistent mapping from CAR RS ids (don't confuse
    with dbSNP rs ids) to more conventional sequence accesssions.

    """

    # MutableMapping doesn't call __missing__ for misses. Emulate that here.
    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            _logger.debug("Cache miss for " + key)
            val = self.__missing__(key)
            super().__setitem__(key, val)
            return val

    def __missing__(self, key):
        url = "http://reg.test.genome.network/refseq/" + key
        resp = requests.get(url)
        if not resp.ok:
            raise KeyError(key)
        ers = resp.json()["externalRecords"]
        for n in ("NCBI", "Ensembl"):
            if n in ers:
                return ers[n]["id"]
        raise KeyError(key)


class AlleleRegistryClient:
    def __init__(self, base_url, login, password, cache_fn=None):
        self._base_url = base_url
        self.login = login
        self.password = password

        if not cache_fn:
            cache_fn = tempfile.NamedTemporaryFile(delete=True).name
        self._refseqmapper = RefSeqMapper(cache_fn)

    def get_allele_by_hgvs(self, hgvs):
        return self._build_vmc_from_ca(self._get("allele", hgvs=hgvs))

    def get_allele(self, hgvs):
        _logger.warn("get_allele is deprecated; use get_allele_by_hgvs instead")
        return self.get_allele_by_hgvs(hgvs)

    def get_allele_by_id(self, id):
        return self._build_vmc_from_ca(self._get("allele/" + id))

    ############################################################################

    def _build_vmc_from_ca(self, cadict):
        def _make_vmc_allele(a):
            """given dict (from CAR json) for single genomicAllele or
            transcriptAllele, create a (Location, Allele) pair, add to
            the bundle, and return the allele.

            """

            car_rsid = a["referenceSequence"].split("/")[-1]
            ir = self._refseqmapper[car_rsid]
            sequence_id = get_vmc_sequence_identifier(ir)

            # N.B. Double check CA coordinate semantics
            # If HGVS like re: insertions, then end -= 1 below
            if len(a["coordinates"]) > 1:
                _logger.warn(f"More than one coordinate set for resp[@id]; using only first")
            coords = a["coordinates"][0]
            interval = models.Interval(start=coords["start"] - 1, end=coords["end"])
            location = models.Location(sequence_id=sequence_id, interval=interval)
            location.id = computed_id(location)
            allele = models.Allele(location_id=location.id, state=coords["allele"])
            allele.id = computed_id(allele)
            return (ir, sequence_id, location, allele)
        
        def _add_ca_to_bundle(bm, ca):
            (ir, sequence_id, location, allele) = _make_vmc_allele(ca)
            bm.locations[location.id] = location
            bm.alleles[allele.id] = allele
            bm.identifiers[sequence_id].update([ir])
            return allele

        # This is the response structure we're going to build
        resp = {
            "@context": None,
            "@id": None,
            "allele_decorations": {},
            "alleles": [],
            "external_records": cadict["externalRecords"],
            "vmcbundle": {},
        }


        alleles = []
        allele_decorations = collections.defaultdict(lambda: dict)
        bm = BundleManager()

        typed_alleles = itertools.chain(
            (("genomic", a) for a in cadict["genomicAlleles"]),
            (("transcript", a) for a in cadict["transcriptAlleles"]))
        for car_type, car_allele in typed_alleles:
            try:
                allele = _add_ca_to_bundle(bm, car_allele)
            except KeyError as e:
                _logger.critical("Failed to make allele for {} ({})".format(
                    str(car_allele), str(e)))
                continue

            allele_decorations[allele.id] = {
                "type": car_type,
                "hgvs": car_allele.get("hgvs", []),
                "gene_id": car_allele.get("geneNCBI_id", None),
                "gene_symbol": car_allele.get("geneSymbol", None),
                "protein_effect": car_allele.get("proteinEffect", {}),
            }
            alleles.append(allele.id)

        resp.update({
            "alleles": alleles,
            "allele_decorations": dict(allele_decorations),
            "vmcbundle": bm.as_dict(),
        })

        return resp

    def _create_credentials(self):
        identity = hashlib.sha1((self.login + self.password).encode('utf-8')).hexdigest()
        gbTime = str(int(time.time()))
        token = hashlib.sha1((url + identity + gbTime).encode('utf-8')).hexdigest()
        return {
            "gbLogin": self.login,
            "gbTime": gbTime,
            "gbToken": token,
        }

    def _get(self, rel_path, **params):
        url = self._base_url + "/" + rel_path
        resp = requests.get(url, params=params)
        resp.raise_for_status()
        return resp.json()



if __name__ == "__main__":
    logging.basicConfig(level="DEBUG")

    config = {
        "base_url": "http://reg.test.genome.network",
        "login": "testuser",
        "password": "testuser"
    }

    arc = AlleleRegistryClient(**config)
    d = arc.get_allele(hgvs="NC_000010.11:g.87894077C>T")
    e = arc._build_vmc_from_ca(d)
    print(json.dumps(e, indent=2, sort_keys=True))
