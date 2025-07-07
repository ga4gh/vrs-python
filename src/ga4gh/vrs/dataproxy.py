"""provides an abstract class for all data access required for
VRS, and a concrete implementation based on seqrepo.

See https://vr-spec.readthedocs.io/en/1.1/impl-guide/required_data.html

"""

import datetime
import functools
import logging
import os
from abc import ABC, abstractmethod
from collections.abc import Sequence
from urllib.parse import urlparse

import requests
from bioutils.accessions import coerce_namespace

_logger = logging.getLogger(__name__)


class DataProxyValidationError(Exception):
    """Class for validation errors during data proxy methods"""


class _DataProxy(ABC):
    """abstract class / interface for VRS data needs

    The proxy MUST support the use of GA4GH sequence identifers (i.e.,
    `ga4gh:SQ...`) as keys, and return these identifiers among the
    aliases for a sequence.  These identifiers may be supported
    natively by the data source or synthesized by the proxy from the
    data source or synthesized.

    """

    @abstractmethod
    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        """Return the specified sequence or subsequence

        start and end are optional

        If the given sequence does not exist, KeyError is raised.

        >> dp.get_sequence("NM_000551.3", 0, 10)
        'CCTCGCCTCC'

        """

    @abstractmethod
    def get_metadata(self, identifier: str) -> dict:
        """For a given identifier, return a structure (dict) containing
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

    @staticmethod
    def extract_sequence_type(alias: str) -> str | None:
        """Provide a convenient way to extract the sequence type from an accession by matching its prefix to a known set of prefixes.

        Args:
        alias (str): The accession string.

        Returns:
        str or None: The sequence type associated with the accession string, or None if no matching prefix is found.

        """
        prefix_dict = {
            "refseq:NM_": "c",
            "refseq:NC_012920": "m",
            "refseq:NG_": "g",
            "refseq:NC_00": "g",
            "refseq:NW_": "g",
            "refseq:NT_": "g",
            "refseq:NR_": "n",
            "refseq:NP_": "p",
            "refseq:XM_": "c",
            "refseq:XR_": "n",
            "refseq:XP_": "p",
            "GRCh": "g",
        }

        for prefix, seq_type in prefix_dict.items():
            if alias.startswith(prefix):
                return seq_type
        return None

    @functools.lru_cache
    def translate_sequence_identifier(
        self, identifier: str, namespace: str | None = None
    ) -> list[str]:
        """Translate given identifier to a list of identifiers in the
        specified namespace.

        `identifier` must be a string
        `namespace` is case-sensitive

        On success, returns string identifier.  Raises KeyError if given
        identifier isn't found.

        """
        try:
            md = self.get_metadata(identifier)
        except (ValueError, KeyError, IndexError) as e:
            raise KeyError(identifier) from e
        aliases = list(set(md["aliases"]))  # ensure uniqueness
        if namespace is not None:
            nsd = namespace + ":"
            aliases = [a for a in aliases if a.startswith(nsd)]
        return aliases

    def derive_refget_accession(self, ac: str) -> str | None:
        """Derive the refget accession from a public accession identifier

        :param ac: public accession in simple or curie form from which to derive the refget accession
        :return: Refget Accession if found
        """
        if ac is None:
            return None

        if ":" not in ac[1:]:
            # always coerce the namespace if none provided
            ac = coerce_namespace(ac)

        refget_accession = None
        try:
            aliases = self.translate_sequence_identifier(ac, namespace="ga4gh")
        except KeyError:
            _logger.exception("KeyError when getting refget accession: %s", ac)
        else:
            if aliases:
                refget_accession = aliases[0].split("ga4gh:")[-1]

        return refget_accession

    def validate_ref_seq(
        self,
        sequence_id: str,
        start_pos: int,
        end_pos: int,
        ref: str,
        require_validation: bool = True,
    ) -> None:
        """Determine whether or not the expected reference sequence matches the actual
        reference sequence. Returns ``None``, but invalid results are logged at level
        WARN by default. If ``require_validation`` is ``True``, then invalid data will
        cause a ``DataProxyValidationError`` to be raised.

        :param sequence_id: Sequence ID to use
        :param start_pos: Start pos (inter-residue) on the sequence_id
        :param end_pos: End pos (inter-residue) on the sequence_id
        :param ref: The expected reference sequence on the sequence_id given the
            start_pos and end_pos
        :param require_validation: If ``True`` and if validation checks fail, a
            ``DataProxyValidationError`` will be raised. Error message will always be
            logged.
        :raises DataProxyValidationError: If excepted reference sequence does not match
            the actual reference sequence and ``require_validation`` is ``True``.
        """
        expected_ref = self.get_sequence(sequence_id, start_pos, end_pos)

        if expected_ref != ref:
            err_msg = f"Reference mismatch at {sequence_id} position {start_pos}-{end_pos} (expected '{expected_ref}' but found '{ref}')"
            _logger.warning(err_msg)

            if require_validation:
                raise DataProxyValidationError(err_msg)


class _SeqRepoDataProxyBase(_DataProxy):
    # wraps seqreqpo classes in order to provide translation to/from
    # `ga4gh` identifiers.

    @functools.lru_cache
    def get_metadata(self, identifier: str) -> dict:
        md = self._get_metadata(identifier)
        md["aliases"] = list(a for a in md["aliases"])  # noqa: C400
        return md

    @functools.lru_cache
    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        return self._get_sequence(identifier, start=start, end=end)

    @abstractmethod
    def _get_metadata(self, identifier: str) -> dict:  # pragma: no cover
        pass

    @abstractmethod
    def _get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:  # pragma: no cover
        pass


class SeqRepoDataProxy(_SeqRepoDataProxyBase):
    """DataProxy based on a local instance of SeqRepo"""

    def __init__(self, sr) -> None:  # noqa: ANN001
        """Initialize DataProxy instance.

        :param sr: SeqRepo instance
        """
        super().__init__()
        self.sr = sr

    def _get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        # fetch raises KeyError if not found
        return self.sr.fetch_uri(coerce_namespace(identifier), start, end)

    def _get_metadata(self, identifier: str) -> dict:
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


class SeqRepoRESTDataProxy(_SeqRepoDataProxyBase):
    """DataProxy based on a REST instance of SeqRepo, as provided by seqrepo-rest-services"""

    rest_version = "1"

    def __init__(self, base_url: str):
        """Initialize REST-based dataproxy instance.

        :param base_url: root URL to server
        """
        super().__init__()
        self.base_url = f"{base_url}/{self.rest_version}/"

    def _get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        url = self.base_url + f"sequence/{identifier}"
        _logger.info("Fetching %s", url)
        params = {"start": start, "end": end}
        resp = requests.get(url, params=params)  # noqa: S113
        if resp.status_code == 404:
            raise KeyError(identifier)
        resp.raise_for_status()
        return resp.text

    def _get_metadata(self, identifier: str) -> dict:
        url = self.base_url + f"metadata/{identifier}"
        _logger.info("Fetching %s", url)
        resp = requests.get(url)  # noqa: S113
        if resp.status_code == 404:
            raise KeyError(identifier)
        resp.raise_for_status()
        return resp.json()


class SequenceProxy(Sequence):
    """Provides efficient and transparent string-like access, including
    random access slicing and reversing, to a biological sequence that
    is stored elsewhere.

    """

    def __init__(self, dp: _DataProxy, alias: str):  # noqa: D107
        self.dp = dp
        self.alias = alias
        self._md = self.dp.get_metadata(self.alias)

    def __str__(self):  # noqa: D105 ANN204
        return self.dp.get_sequence(self.alias)

    def __len__(self):  # noqa: D105 ANN204
        return self._md["length"]

    def __reversed__(self):  # noqa: D105 ANN204
        msg = "Reversed iteration of a SequenceProxy is not implemented"
        raise NotImplementedError(msg)

    def __getitem__(self, key):  # noqa: ANN001 ANN204
        """Return sequence for key (slice), fetching if necessary"""
        if isinstance(key, int):
            key = slice(key, key + 1)
        if key.step is not None:
            msg = "Only contiguous sequence slices are supported"
            raise ValueError(msg)

        return self.dp.get_sequence(self.alias, key.start, key.stop)


def _isoformat(o: datetime.datetime) -> str:
    """Convert datetime.datetime to iso formatted timestamp

    >>> dt = datetime.datetime(2019, 10, 15, 10, 23, 41, 115927)
    >>> _isoformat(dt)
    '2019-10-15T10:23:41.115927Z'

    :param o: date object to format
    :return: ISO8601-formatted string equivalent
    :raise TypeError: if object given isn't a Python datetime instance
    """
    # stolen from connexion flask_app.py
    if not isinstance(o, datetime.datetime):
        raise TypeError
    if o.tzinfo:
        # eg: '2015-09-25T23:14:42.588601+00:00'
        return o.isoformat("T")
    # No timezone present - assume UTC.
    # eg: '2015-09-25T23:14:42.588601Z'
    return o.isoformat("T") + "Z"


# Future implementations
# * The RefGetDataProxy is waiting on support for sequence lookup by alias
# class RefGetDataProxy(_DataProxy):
#     def __init__(self, base_url):
#         super().__init__()
#         self.base_url = base_url


def create_dataproxy(uri: str | None = None) -> _DataProxy:
    """Create a dataproxy from uri or GA4GH_VRS_DATAPROXY_URI

    Currently accepted URI schemes:

    * seqrepo+file:///path/to/seqrepo/root
    * seqrepo+:../relative/path/to/seqrepo/root
    * seqrepo+http://localhost:5000/seqrepo
    * seqrepo+https://somewhere:5000/seqrepo

    :raise ValueError: if URI doesn't match recognized schemes, e.g. is missing provider
        prefix (`"seqrepo+"`)
    """
    uri = uri or os.environ.get("GA4GH_VRS_DATAPROXY_URI", None)

    if uri is None:
        msg = "No data proxy URI provided or found in GA4GH_VRS_DATAPROXY_URI"
        raise ValueError(msg)

    parsed_uri = urlparse(uri)
    scheme = parsed_uri.scheme

    if "+" not in scheme:
        msg = "create_dataproxy scheme must include provider (e.g., `seqrepo+http:...`)"
        raise ValueError(msg)

    provider, proto = scheme.split("+")

    if provider == "seqrepo":
        if proto in ("", "file"):
            from biocommons.seqrepo import SeqRepo

            sr = SeqRepo(root_dir=parsed_uri.path)
            dp = SeqRepoDataProxy(sr)
        elif proto in ("http", "https"):
            dp = SeqRepoRESTDataProxy(uri[len(provider) + 1 :])
        else:
            msg = f"SeqRepo URI scheme {parsed_uri.scheme} not implemented"
            raise ValueError(msg)

    else:
        msg = f"DataProxy provider {provider} not implemented"
        raise ValueError(msg)

    return dp
