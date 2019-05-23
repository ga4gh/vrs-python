"""provides an abstract class for all data access required for
vr.extras, and a concrete implementation based on seqrepo.

"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime
from typing import List

from ga4gh.vr.utils import coerce_namespace


@dataclass
class SequenceMetadata:
    # intentionally hidden: sequence_id: str
    length: int
    alphabet: str
    added: datetime
    aliases: List[str]


class VRDataProxy(ABC):
    @abstractmethod
    def get_sequence(identifier, start=None, end=None):
        """return the specified sequence or subsequence

        start and end are optional

        """

    @abstractmethod
    def get_metadata(identifier):
        """for a given identifier, return a structure (dict) containing
        sequence length, aliases, and other optional info

        """


class SeqRepoDataProxy(VRDataProxy):
    def __init__(self, sr):
        super().__init__()
        self.sr = sr

    def get_sequence(self, alias, start=None, end=None):
        return self.sr.fetch(alias, start, end)

    def get_metadata(self, alias):
        ns, a = coerce_namespace(alias).split(":", 2)
        ns = "RefSeq" if ns == "refseq" else ns
        r = self.sr.aliases.find_aliases(namespace=ns, alias=a).fetchone()
        seqinfo = self.sr.sequences.fetch_seqinfo(r["seq_id"])
        aliases = self.sr.aliases.fetch_aliases(r["seq_id"])
        md = SequenceMetadata(
            length=seqinfo["len"],
            alphabet=seqinfo["alpha"],
            added=seqinfo["added"],
            aliases=[f"{a['namespace']}:{a['alias']}" for a in aliases],
            )
        return md


# Future implementations
# * The SeqRepoRESTDataProxy is waiting on the SeqRepo REST interface
# * The RefGetDataProxy is waiting on support for sequence lookup by alias
# class SeqRepoRESTDataProxy(VRDataProxy):
#     def __init__(self, base_url):
#         super().__init__()
#         self.base_url = base_url
#  
#  
# class RefGetDataProxy(VRDataProxy):
#     def __init__(self, base_url):
#         super().__init__()
#         self.base_url = base_url

