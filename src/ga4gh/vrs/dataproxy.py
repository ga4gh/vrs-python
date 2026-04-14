"""Provides abstract interfaces for all data access required by VRS,
and a concrete implementation based on seqrepo.

See https://vr-spec.readthedocs.io/en/1.1/impl-guide/required_data.html

"""

from __future__ import annotations

import functools
import logging
import os
from abc import ABC, abstractmethod
from collections.abc import Callable, Sequence
from urllib.parse import urlparse

_logger = logging.getLogger(__name__)


class DataProxyValidationError(Exception):
    """Class for validation errors during data proxy methods"""


class SequenceDataProxy(ABC):
    """Abstract interface for sequence data access.

    Provides methods for retrieving genomic sequences and their lengths.
    """

    @abstractmethod
    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        """Return the specified sequence or subsequence.

        start and end are optional.

        If the given sequence does not exist, KeyError is raised.

        >> dp.get_sequence("NM_000551.3", 0, 10)
        'CCTCGCCTCC'
        """

    @abstractmethod
    def get_sequence_length(self, identifier: str) -> int:
        """Return the length of the specified sequence.

        If the given sequence does not exist, KeyError is raised.

        >> dp.get_sequence_length("NM_000551.3")
        4560
        """

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
        correct_ref = self.get_sequence(sequence_id, start_pos, end_pos)

        if correct_ref != ref:
            err_msg = f"Reference mismatch at {sequence_id} position {start_pos}-{end_pos} (input gave '{ref}' but correct ref is '{correct_ref}')"
            _logger.warning(err_msg)

            if require_validation:
                raise DataProxyValidationError(err_msg)


class AliasDataProxy(ABC):
    """Abstract interface for sequence alias / identifier resolution.

    Provides methods for translating between different sequence identifier
    namespaces and deriving refget accessions.
    """

    @abstractmethod
    def get_aliases(self, identifier: str) -> list[str]:
        """Return a list of aliases for the given sequence identifier.

        Each alias is a string of the form ``namespace:accession``.

        If the given sequence does not exist, KeyError is raised.

        >> dp.get_aliases("NM_000551.3")
        ['MD5:215137b1973c1a5afcf86be7d999574a',
         'RefSeq:NM_000551.3',
         'SEGUID:T12L0p2X5E8DbnL0+SwI4Wc1S6g',
         'SHA1:4f5d8bd29d97e44f036e72f4f92c08e167354ba8',
         'ga4gh:SQ.v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_',
         'gi:319655736']
        """

    def translate_sequence_identifier(
        self, identifier: str, namespace: str | None = None
    ) -> list[str]:
        """Translate given identifier to a list of identifiers in the
        specified namespace.

        ``identifier`` must be a string.
        ``namespace`` is case-sensitive.

        On success, returns list of string identifiers. Raises KeyError if given
        identifier isn't found.
        """
        try:
            aliases = self.get_aliases(identifier)
        except (ValueError, KeyError, IndexError) as e:
            raise KeyError(identifier) from e
        aliases = list(set(aliases))  # ensure uniqueness
        if namespace is not None:
            nsd = namespace + ":"
            aliases = [a for a in aliases if a.startswith(nsd)]
        return aliases

    def derive_refget_accession(self, ac: str) -> str | None:
        """Derive the refget accession from a public accession identifier.

        :param ac: public accession in simple or curie form from which to derive the refget accession
        :return: Refget Accession if found
        """
        if ac is None:
            return None

        if ":" not in ac[1:]:
            # always coerce the namespace if none provided
            from bioutils.accessions import coerce_namespace

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


class DataProxy(SequenceDataProxy, AliasDataProxy):
    """Combined abstract interface for sequence data access and alias mapping.

    The proxy MUST support the use of GA4GH sequence identifiers (i.e.,
    ``ga4gh:SQ...``) as keys, and return these identifiers among the
    aliases for a sequence. These identifiers may be supported
    natively by the data source or synthesized by the proxy from the
    data source.
    """


class CompositeDataProxy(DataProxy):
    """Composes a SequenceDataProxy and AliasDataProxy into a full DataProxy.

    For sequence methods: tries the sequence backend first; on KeyError,
    translates the identifier to a GA4GH digest via the alias backend and retries.
    For alias methods: delegates directly to the alias backend.
    """

    def __init__(
        self, sequence_proxy: SequenceDataProxy, alias_proxy: AliasDataProxy
    ) -> None:
        """Initialize CompositeDataProxy.

        :param sequence_proxy: Backend for sequence retrieval
        :param alias_proxy: Backend for alias / identifier resolution
        """
        super().__init__()
        self._sequence_proxy = sequence_proxy
        self._alias_proxy = alias_proxy

    @functools.lru_cache(maxsize=512)
    def _resolve_ga4gh_identifier(self, identifier: str) -> str:
        """Translate identifier to a GA4GH SQ digest via the alias backend.

        :param identifier: Sequence identifier in any namespace
        :return: GA4GH identifier (``ga4gh:SQ.…``)
        :raises KeyError: if no GA4GH alias is found
        """
        aliases = self._alias_proxy.translate_sequence_identifier(
            identifier, namespace="ga4gh"
        )
        sq_aliases = [a for a in aliases if a.startswith("ga4gh:SQ.")]
        if not sq_aliases:
            raise KeyError(identifier)
        return sq_aliases[0]

    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        """Return the specified sequence or subsequence.

        Tries the sequence backend first; on KeyError, resolves the identifier
        to a GA4GH digest via the alias backend and retries.
        """
        try:
            return self._sequence_proxy.get_sequence(identifier, start, end)
        except KeyError:
            pass
        ga4gh_id = self._resolve_ga4gh_identifier(identifier)
        return self._sequence_proxy.get_sequence(ga4gh_id, start, end)

    def get_sequence_length(self, identifier: str) -> int:
        """Return the length of the specified sequence.

        Tries the sequence backend first; on KeyError, resolves the identifier
        to a GA4GH digest via the alias backend and retries.
        """
        try:
            return self._sequence_proxy.get_sequence_length(identifier)
        except KeyError:
            pass
        ga4gh_id = self._resolve_ga4gh_identifier(identifier)
        return self._sequence_proxy.get_sequence_length(ga4gh_id)

    def get_aliases(self, identifier: str) -> list[str]:
        """Return all aliases for the specified sequence, delegating to the alias backend."""
        return self._alias_proxy.get_aliases(identifier)


class SequenceProxy(Sequence):
    """Provides efficient and transparent string-like access, including
    random access slicing and reversing, to a biological sequence that
    is stored elsewhere.

    """

    def __init__(self, dp: SequenceDataProxy, alias: str) -> None:  # noqa: D107
        self.dp = dp
        self.alias = alias
        self._length = self.dp.get_sequence_length(self.alias)

    def __str__(self) -> str:  # noqa: D105
        return self.dp.get_sequence(self.alias)

    def __len__(self):  # noqa: D105 ANN204
        return self._length

    def __reversed__(self):  # noqa: D105 ANN204
        msg = "Reversed iteration of a SequenceProxy is not implemented"
        raise NotImplementedError(msg)

    def __getitem__(self, key) -> str:  # noqa: ANN001
        """Return sequence for key (slice), fetching if necessary"""
        if isinstance(key, int):
            key = slice(key, key + 1)
        if key.step is not None:
            msg = "Only contiguous sequence slices are supported"
            raise ValueError(msg)

        return self.dp.get_sequence(self.alias, key.start, key.stop)


# Provider registry for extensible create_dataproxy
_provider_registry: dict[str, Callable[..., DataProxy]] = {}


def register_provider(name: str, factory: Callable[..., DataProxy]) -> None:
    """Register a custom DataProxy provider factory.

    :param name: Provider name used in URIs (e.g. "seqrepo" for "seqrepo+http://...")
    :param factory: Callable that accepts (parsed_uri, proto, disable_healthcheck) and returns a DataProxy
    """
    _provider_registry[name] = factory


def create_dataproxy(
    uri: str | None = None, disable_healthcheck: bool = False
) -> DataProxy:
    """Create a dataproxy from uri or GA4GH_VRS_DATAPROXY_URI

    :param uri: Dataproxy URI. Currently accepted URI schemes:

        * seqrepo+file:///path/to/seqrepo/root
        * seqrepo+:../relative/path/to/seqrepo/root
        * seqrepo+http://localhost:5000/seqrepo
        * seqrepo+https://somewhere:5000/seqrepo
        * refget+file:///path/to/refgetstore/root
    :param disable_healthcheck: Whether healthcheck should be disabled in REST dataproxy
    :raise ValueError: if URI doesn't match recognized schemes, e.g. is missing provider
        prefix (``"seqrepo+"``)
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

    if provider in _provider_registry:
        return _provider_registry[provider](
            parsed_uri, proto, disable_healthcheck=disable_healthcheck
        )

    if provider == "seqrepo":
        from ga4gh.vrs.dataproxy_seqrepo import _create_seqrepo_dataproxy

        return _create_seqrepo_dataproxy(
            uri, parsed_uri, proto, disable_healthcheck=disable_healthcheck
        )

    if provider == "refget":
        from ga4gh.vrs.dataproxy_refget import _create_refget_dataproxy

        return _create_refget_dataproxy(
            parsed_uri, proto, disable_healthcheck=disable_healthcheck
        )

    msg = f"DataProxy provider {provider} not implemented"
    raise ValueError(msg)
