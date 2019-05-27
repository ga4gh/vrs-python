"""provides an abstract class for all data access required for
vr.extras, and a concrete implementation based on seqrepo.

"""

from abc import ABC, abstractmethod
from datetime import datetime
import logging

from bioutils.accessions import coerce_namespace
import requests

from .utils import isoformat


_logger = logging.getLogger(__name__)


class _DataProxy(ABC):
    """abstract class / interface for vr data needs
    
    """

    @abstractmethod
    def get_sequence(identifier, start=None, end=None):
        """return the specified sequence or subsequence

        start and end are optional

        If the given sequence does not exist, KeyError is raised.

        """

    @abstractmethod
    def get_metadata(identifier):
        """for a given identifier, return a structure (dict) containing
        sequence length, aliases, and other optional info

        If the given sequence does not exist, KeyError is raised.

        """


class SeqRepoDataProxy(_DataProxy):
    def __init__(self, sr):
        super().__init__()
        self.sr = sr

    def get_sequence(self, identifier, start=None, end=None):
        # fetch raises KeyError if not found
        return self.sr.fetch(identifier, start, end)
        

    def get_metadata(self, identifier):
        ns, a = coerce_namespace(identifier).split(":", 2)
        ns = "RefSeq" if ns == "refseq" else ns
        r = self.sr.aliases.find_aliases(namespace=ns, alias=a).fetchone()
        if r is None:
            raise KeyError(identifier) 
        seqinfo = self.sr.sequences.fetch_seqinfo(r["seq_id"])
        aliases = self.sr.aliases.fetch_aliases(r["seq_id"])
        return {
            "length": seqinfo["len"],
            "alphabet": seqinfo["alpha"],
            "added": isoformat(seqinfo["added"]),
            "aliases": [f"{a['namespace']}:{a['alias']}" for a in aliases],
            }
        return md


class SeqRepoRESTDataProxy(_DataProxy):
    rest_version = "1"

    def __init__(self, base_url):
        super().__init__()
        self.base_url = f"{base_url}/{self.rest_version}/"

    def get_sequence(self, identifier, start=None, end=None):
        url = self.base_url + f"sequence/{identifier}"
        _logger.info("Fetching " + url)
        params = {"start": start, "end": end}
        resp = requests.get(url, params=params)
        if resp.status_code == 404:
            raise KeyError(identifier) 
        resp.raise_for_status()
        return resp.text

    def get_metadata(self, identifier):
        url = self.base_url + f"metadata/{identifier}"
        _logger.info("Fetching " + url)
        resp = requests.get(url)
        if resp.status_code == 404:
            raise KeyError(identifier) 
        resp.raise_for_status()
        data = resp.json()["metadata"]
        return data
    

 
# Future implementations
# * The RefGetDataProxy is waiting on support for sequence lookup by alias
# class RefGetDataProxy(_DataProxy):
#     def __init__(self, base_url):
#         super().__init__()
#         self.base_url = base_url

