"""provides an abstract class for all data access required for
VRS, and a concrete implementation based on seqrepo.

See https://vr-spec.readthedocs.io/en/1.1/impl-guide/required_data.html

"""

from abc import ABC, abstractmethod
from collections.abc import Sequence
import datetime
import functools
import logging
import os
from urllib.parse import urlparse

from bioutils.accessions import coerce_namespace
import requests



_logger = logging.getLogger(__name__)


class _DataProxy(ABC):
    """abstract class / interface for VRS data needs
    
    The proxy MUST support the use of GA4GH sequence identifers (i.e.,
    `ga4gh:SQ...`) as keys, and return these identifiers among the
    aliases for a sequence.  These identifiers may be supported
    natively by the data source or synthesized by the proxy from the
    data source or synthesized.

    """

    @abstractmethod
    def get_sequence(self, identifier, start=None, end=None):
        """return the specified sequence or subsequence

n        start and end are optional

        If the given sequence does not exist, KeyError is raised.

        >> dp.get_sequence("NM_000551.3", 0, 10)
        'CCTCGCCTCC'

        """

    @abstractmethod
    def get_metadata(self, identifier):
        """for a given identifier, return a structure (dict) containing
        sequence length, aliases, and other optional info

        If the given sequence does not exist, KeyError is raised.

        >> dp.get_metadata("NM_000551.3")
        {'added': '2016-08-24T05:03:11Z',
         'aliases': ['MD5:215137b1973c1a5afcf86be7d999574a',
                     'RefSeq:NM_000551.3',
                     'SEGUID:T12L0p2X5E8DbnL0+SwI4Wc1S6g',
                     'SHA1:4f5d8bd29d97e44f036e72f4f92c08e167354ba8',
                     'ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_',
                     'gi:319655736'],
         'alphabet': 'ACGT',
         'length': 4560}

        """


    @functools.lru_cache()
    def translate_sequence_identifier(self, identifier, namespace=None):
        """Translate given identifier to a list of identifiers in the
        specified namespace.

        `identifier` must be a string
        `namespace` is case-sensitive

        On success, returns string identifier.  Raises KeyError if given
        identifier isn't found.

        """

        try:
            md = self.get_metadata(identifier)
        except (ValueError, KeyError, IndexError):
            raise KeyError(identifier)
        aliases = list(set(md["aliases"]))  # ensure uniqueness
        if namespace is not None:
            nsd = namespace + ":"
            aliases = [a for a in aliases if a.startswith(nsd)]
        return aliases


class _SeqRepoDataProxyBase(_DataProxy):
    # wraps seqreqpo classes in order to provide translation to/from
    # `ga4gh` identifiers.

    def get_metadata(self, identifier):
        md = self._get_metadata(identifier)
        md["aliases"] = list(a for a in md["aliases"])
        return md

    def get_sequence(self, identifier, start=None, end=None):
        return self._get_sequence(identifier, start=start, end=end)

    @abstractmethod
    def _get_metadata(self, identifier):  # pragma: no cover
        pass

    @abstractmethod
    def _get_sequence(self, identifier, start=None, end=None):  # pragma: no cover
        pass


class SeqRepoDataProxy(_SeqRepoDataProxyBase):
    """DataProxy based on a local instance of SeqRepo"""

    def __init__(self, sr):
        super().__init__()
        self.sr = sr

    def _get_sequence(self, identifier, start=None, end=None):
        # fetch raises KeyError if not found
        return self.sr.fetch_uri(identifier, start, end)

    def _get_metadata(self, identifier):
        ns, a = coerce_namespace(identifier).split(":", 2)
        r = list(self.sr.aliases.find_aliases(namespace=ns, alias=a))
        if len(r) == 0:
            raise KeyError(identifier) 
        seq_id = r[0]["seq_id"]
        seqinfo = self.sr.sequences.fetch_seqinfo(seq_id)
        aliases = self.sr.aliases.find_aliases(seq_id=seq_id)
        return {
            "length": seqinfo["len"],
            "alphabet": seqinfo["alpha"],
            "added": _isoformat(seqinfo["added"]),
            "aliases": [f"{a['namespace']}:{a['alias']}" for a in aliases],
            }
        return md


class SeqRepoRESTDataProxy(_SeqRepoDataProxyBase):
    """DataProxy based on a REST instance of SeqRepo, as provided by seqrepo-rest-services"""

    rest_version = "1"

    def __init__(self, base_url):
        super().__init__()
        self.base_url = f"{base_url}/{self.rest_version}/"

    def _get_sequence(self, identifier, start=None, end=None):
        url = self.base_url + f"sequence/{identifier}"
        _logger.info("Fetching " + url)
        params = {"start": start, "end": end}
        resp = requests.get(url, params=params)
        if resp.status_code == 404:
            raise KeyError(identifier) 
        resp.raise_for_status()
        return resp.text

    def _get_metadata(self, identifier):
        url = self.base_url + f"metadata/{identifier}"
        _logger.info("Fetching " + url)
        resp = requests.get(url)
        if resp.status_code == 404:
            raise KeyError(identifier) 
        resp.raise_for_status()
        data = resp.json()
        return data
    

class SequenceProxy(Sequence):
    """Provides efficient and transparent string-like access, including
    random access slicing and reversing, to a biological sequence that
    is stored elsewhere.

    """

    def __init__(self, dp, alias):
        self.dp = dp
        self.alias = alias
        self._md = self.dp.get_metadata(self.alias)
        
    def __str__(self):
        return self.dp.get_sequence(self.alias)

    def __len__(self):
        return self._md["length"]

    def __reversed__(self):
        raise NotImplementedError("Reversed iteration of a SequenceProxy is not implemented")

    def __getitem__(self, key):
        """return sequence for key (slice), fetching if necessary

        """

        if isinstance(key, int):
            key = slice(key, key+1)
        if key.step is not None:
            raise ValueError("Only contiguous sequence slices are supported")

        return self.dp.get_sequence(self.alias, key.start, key.stop)



def _isoformat(o):
    """convert datetime.datetime to iso formatted timestamp

    >>> dt = datetime.datetime(2019, 10, 15, 10, 23, 41, 115927)
    >>> _isoformat(dt)
    '2019-10-15T10:23:41.115927Z'

    """

    # stolen from connexion flask_app.py
    assert isinstance(o, datetime.datetime)
    if o.tzinfo:
        # eg: '2015-09-25T23:14:42.588601+00:00'
        return o.isoformat('T')
    # No timezone present - assume UTC.
    # eg: '2015-09-25T23:14:42.588601Z'
    return o.isoformat('T') + 'Z'

 
# Future implementations
# * The RefGetDataProxy is waiting on support for sequence lookup by alias
# class RefGetDataProxy(_DataProxy):
#     def __init__(self, base_url):
#         super().__init__()
#         self.base_url = base_url



def create_dataproxy(uri: str = None) -> _DataProxy:
    """Create a dataproxy from uri or GA4GH_VRS_DATAPROXY_URI
    
    Currently accepted URI schemes:

    * seqrepo+file:///path/to/seqrepo/root
    * seqrepo+:../relative/path/to/seqrepo/root
    * seqrepo+http://localhost:5000/seqrepo
    * seqrepo+https://somewhere:5000/seqrepo

    """

    uri = (uri
           or os.environ.get("GA4GH_VRS_DATAPROXY_URI", None))

    if uri is None:
        raise ValueError("No data proxy URI provided or found in GA4GH_VRS_DATAPROXY_URI")

    parsed_uri = urlparse(uri)
    scheme = parsed_uri.scheme
    
    if "+" not in scheme:
        raise ValueError("create_dataproxy scheme must include provider (e.g., `seqrepo+http:...`)")

    provider, proto = scheme.split("+")

    if provider == "seqrepo":
        if proto in ("", "file"):
            from biocommons.seqrepo import SeqRepo
            sr = SeqRepo(root_dir=parsed_uri.path)
            dp = SeqRepoDataProxy(sr)
        elif proto in ("http", "https"):
            dp = SeqRepoRESTDataProxy(uri[len(provider)+1:])
        else:
            raise ValueError(f"SeqRepo URI scheme {parsed_uri.scheme} not implemented")

    else:
        raise ValueError(f"DataProxy provider {provider} not implemented")

    return dp




if __name__ == "__main__":
    # Before running, do something like this:
    # snafu$ docker run \
    # >   --name seqrepo-rest-service \
    # >   --detach --rm -p 5000:5000 \
    # >   -v /usr/local/share/seqrepo/:/usr/local/share/seqrepo/ \
    # >   biocommons/seqrepo-rest-service

    dp1 = create_dataproxy("seqrepo+http://localhost:5000/seqrepo")
    dp2 = create_dataproxy("seqrepo+file:///usr/local/share/seqrepo/latest")

    ir = "refseq:NM_000551.3"

    print(f"dp1 = {dp1}")
    print(f"dp2 = {dp2}")

    assert dp1.get_metadata(ir) == dp2.get_metadata(ir)
    assert dp1.get_sequence(ir) == dp2.get_sequence(ir)
    
