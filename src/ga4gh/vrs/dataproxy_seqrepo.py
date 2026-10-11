"""SeqRepo-based implementations of the DataProxy interface.

Provides concrete DataProxy implementations backed by a local SeqRepo
instance or a SeqRepo REST service.
"""

from __future__ import annotations

import datetime
import functools
import logging
from urllib.parse import ParseResult

import requests
from bioutils.accessions import coerce_namespace

from ga4gh.vrs.dataproxy import DataProxy

_logger = logging.getLogger(__name__)


class SeqRepoDataProxy(DataProxy):
    """DataProxy based on a local instance of SeqRepo."""

    def __init__(self, sr) -> None:  # noqa: ANN001
        """Initialize DataProxy instance.

        :param sr: SeqRepo instance
        """
        super().__init__()
        self.sr = sr

    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        """Return the specified sequence or subsequence from local SeqRepo."""
        return self.sr.fetch_uri(coerce_namespace(identifier), start, end)

    @functools.lru_cache(maxsize=512)  # noqa: B019
    def get_sequence_length(self, identifier: str) -> int:
        """Return the length of the specified sequence from local SeqRepo."""
        ns, a = coerce_namespace(identifier).split(":", 2)
        r = list(self.sr.aliases.find_aliases(namespace=ns, alias=a))
        if len(r) == 0:
            raise KeyError(identifier)
        seq_id = r[0]["seq_id"]
        seqinfo = self.sr.sequences.fetch_seqinfo(seq_id)
        return seqinfo["len"]

    @functools.lru_cache(maxsize=512)  # noqa: B019
    def get_aliases(self, identifier: str) -> list[str]:
        """Return all aliases for the specified sequence from local SeqRepo."""
        ns, a = coerce_namespace(identifier).split(":", 2)
        r = list(self.sr.aliases.find_aliases(namespace=ns, alias=a))
        if len(r) == 0:
            raise KeyError(identifier)
        seq_id = r[0]["seq_id"]
        aliases = self.sr.aliases.find_aliases(seq_id=seq_id)
        return [f"{a['namespace']}:{a['alias']}" for a in aliases]


class SeqRepoRESTDataProxy(DataProxy):
    """DataProxy based on a REST instance of SeqRepo, as provided by seqrepo-rest-services."""

    rest_version = "1"

    def __init__(self, base_url: str, disable_healthcheck: bool = False) -> None:
        """Initialize REST-based dataproxy instance.

        :param base_url: root URL to server
        :param disable_healthcheck: Whether healthcheck should be disabled
        """
        super().__init__()
        self.base_url = f"{base_url}/{self.rest_version}/"
        if not disable_healthcheck:
            ping_url = self.base_url + "ping"
            ping_resp = requests.get(ping_url)  # noqa: S113
            ping_resp.raise_for_status()

    @functools.lru_cache  # noqa: B019
    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        """Return the specified sequence or subsequence from REST SeqRepo."""
        url = self.base_url + f"sequence/{identifier}"
        _logger.info("Fetching %s", url)
        params = {"start": start, "end": end}
        resp = requests.get(url, params=params)  # noqa: S113
        if resp.status_code == 404:
            raise KeyError(identifier)
        resp.raise_for_status()
        return resp.text

    @functools.lru_cache  # noqa: B019
    def _get_metadata(self, identifier: str) -> dict:
        """Fetch metadata from the REST service.

        Internal helper that returns the full metadata dict from the REST API.
        Used by both get_sequence_length and get_aliases to avoid duplicate requests.
        """
        url = self.base_url + f"metadata/{identifier}"
        _logger.info("Fetching %s", url)
        resp = requests.get(url)  # noqa: S113
        if resp.status_code == 404:
            raise KeyError(identifier)
        resp.raise_for_status()
        return resp.json()

    def get_sequence_length(self, identifier: str) -> int:
        """Return the length of the specified sequence from REST SeqRepo."""
        md = self._get_metadata(identifier)
        return md["length"]

    def get_aliases(self, identifier: str) -> list[str]:
        """Return all aliases for the specified sequence from REST SeqRepo."""
        md = self._get_metadata(identifier)
        return list(md["aliases"])

    def get_metadata(self, identifier: str) -> dict:
        """Return length and aliases in a single REST round-trip.

        Overrides the default :meth:`DataProxy.get_metadata` to avoid the two
        sequential calls that would otherwise each hit ``_get_metadata``.
        """
        md = self._get_metadata(identifier)
        return {"length": md["length"], "aliases": list(md["aliases"])}


def _isoformat(o: datetime.datetime) -> str:
    """Convert datetime.datetime to iso formatted timestamp.

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


def _create_seqrepo_dataproxy(
    uri: str,
    parsed_uri: ParseResult,
    proto: str,
    *,
    disable_healthcheck: bool = False,
) -> DataProxy:
    """Create a SeqRepo-backed DataProxy from a parsed URI.

    :param uri: Original URI string
    :param parsed_uri: Parsed URI components
    :param proto: Protocol portion of the scheme (e.g. "file", "http", "https", "")
    :param disable_healthcheck: Whether healthcheck should be disabled in REST dataproxy
    :return: A SeqRepoDataProxy or SeqRepoRESTDataProxy instance
    """
    if proto in ("", "file"):
        from biocommons.seqrepo import SeqRepo

        sr = SeqRepo(root_dir=parsed_uri.path)
        return SeqRepoDataProxy(sr)

    if proto in ("http", "https"):
        provider = parsed_uri.scheme.split("+")[0]
        return SeqRepoRESTDataProxy(
            uri[len(provider) + 1 :], disable_healthcheck=disable_healthcheck
        )

    msg = f"SeqRepo URI scheme {parsed_uri.scheme} not implemented"
    raise ValueError(msg)
