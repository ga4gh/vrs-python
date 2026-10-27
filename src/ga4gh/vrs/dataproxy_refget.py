"""Refget/gtars-based implementations of the DataProxy interfaces.

Provides a ``SequenceDataProxy`` backed by a gtars ``RefgetStore`` for
high-performance sequence retrieval by GA4GH digest, and a matching
``AliasDataProxy`` that uses the store's sequence alias index to resolve
cross-namespace identifiers. The module factory composes both into a full
``DataProxy`` from a single URI. Stores may be opened from a local directory
or from a remote static HTTP(S) tree. Remote stores use a local cache rooted
at ``GA4GH_VRS_REFGET_CACHE_DIR`` when set, otherwise under an XDG-style
default cache directory.

Requires the ``refget`` optional dependency group: ``pip install ga4gh.vrs[refget]``

Known divergence from SeqRepo
-----------------------------
For protein sequences containing embedded stop-codon ``*`` characters, a
RefgetStore can yield a different ``ga4gh:SQ.`` digest than SeqRepo for the
*same underlying bytes*. SeqRepo runs the FASTA through
``bioutils.sequences.normalize_sequence()``, which strips ``*`` before
computing the digest; gtars hashes the raw bytes. The RefgetStore digest is
the spec-correct one (VRS does not prescribe ``*`` stripping). In practice
this only affects a small fraction (~0.1%) of Ensembl protein accessions;
NCBI RefSeq protein FASTAs contain no ``*`` and DNA/RNA sequences are never
affected. See ``docs/dataproxy_backends.md`` for details.
"""

from __future__ import annotations

import functools
import hashlib
import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING
from urllib.parse import ParseResult, urlunparse

from ga4gh.vrs.dataproxy import (
    DATAPROXY_CACHE_SIZE,
    AliasDataProxy,
    CompositeDataProxy,
    DataProxy,
    SequenceDataProxy,
    coerce_accession_namespace,
    validate_start_le_end,
)

if TYPE_CHECKING:
    from gtars.refget import RefgetStore

_logger = logging.getLogger(__name__)

_GA4GH_SQ_PREFIX = "ga4gh:SQ."
_SQ_PREFIX = "SQ."
_SHA512T24U_PREFIX = "sha512t24u:"
_REFGET_CACHE_ENV_VAR = "GA4GH_VRS_REFGET_CACHE_DIR"


def _extract_digest(identifier: str) -> str | None:
    """Strip GA4GH SQ prefix, returning the raw sha512t24u digest.

    Accepts identifiers in ``ga4gh:SQ.<digest>``, ``SQ.<digest>``, or
    ``sha512t24u:<digest>`` form.

    :param identifier: Sequence identifier, possibly with GA4GH prefix
    :return: Raw digest string, or ``None`` if identifier is not in GA4GH SQ form
    """
    if identifier.startswith(_GA4GH_SQ_PREFIX):
        return identifier[len(_GA4GH_SQ_PREFIX) :]
    if identifier.startswith(_SQ_PREFIX):
        return identifier[len(_SQ_PREFIX) :]
    if identifier.startswith(_SHA512T24U_PREFIX):
        return identifier[len(_SHA512T24U_PREFIX) :]
    return None


def _slug_remote_url(remote_url: str) -> str:
    """Return a short deterministic cache slug for a remote refget URL."""
    return hashlib.sha256(remote_url.encode("utf-8")).hexdigest()[:16]


def _resolve_refget_cache_base_dir() -> Path:
    """Return the base directory for remote refget caches.

    This applies only to remote ``refget+http://...`` and
    ``refget+https://...`` URIs. Resolution order is:

    1. ``GA4GH_VRS_REFGET_CACHE_DIR``, if set
    2. ``$XDG_CACHE_HOME/gtars/refget``, if ``XDG_CACHE_HOME`` is set
    3. ``~/.cache/gtars/refget``

    The returned path is a cache *root*. The actual cache directory passed to
    ``RefgetStore.open_remote()`` is a deterministic per-remote subdirectory
    under this root; see :func:`_resolve_remote_cache_dir`.
    """
    cache_root = os.environ.get(_REFGET_CACHE_ENV_VAR)
    if cache_root:
        return Path(cache_root).expanduser()

    xdg_cache_home = os.environ.get("XDG_CACHE_HOME")
    if xdg_cache_home:
        return Path(xdg_cache_home).expanduser() / "gtars" / "refget"

    return Path.home() / ".cache" / "gtars" / "refget"


def _resolve_remote_cache_dir(remote_url: str) -> Path:
    """Return the per-remote cache directory for a remote refget URL."""
    return _resolve_refget_cache_base_dir() / _slug_remote_url(remote_url)


class RefgetSequenceDataProxy(SequenceDataProxy):
    """SequenceDataProxy backed by a gtars RefgetStore.

    Only accepts ``ga4gh:SQ.<digest>`` or ``SQ.<digest>`` identifiers.
    Non-GA4GH identifiers raise KeyError — use CompositeDataProxy with an
    AliasDataProxy to resolve namespace identifiers first.
    """

    def __init__(self, store: RefgetStore) -> None:
        """Initialize RefgetSequenceDataProxy.

        :param store: A gtars RefgetStore instance (local, on-disk, or in-memory)
        """
        super().__init__()
        self.store = store
        # Track which digests we've already promoted from Stub to Full so
        # subsequent substring fetches skip the per-digest file read.
        self._loaded_digests: set[str] = set()

    def _ensure_sequence_loaded(self, digest: str) -> None:
        """Lazily materialize sequence bytes for ``digest`` in the store.

        ``open_local()`` reads ``sequences.rgsi`` and inserts a Stub record
        for every sequence in the store, so metadata lookups work without
        any preload. But byte-level access (``get_substring``) needs the
        per-digest payload read in via ``load_sequence(digest)``. We call
        that on first access for each digest and memoize the result — no
        sequence bytes are loaded until a caller actually asks for them.
        """
        if digest in self._loaded_digests:
            return
        self.store.load_sequence(digest)  # type: ignore[attr-defined]
        self._loaded_digests.add(digest)

    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        """Return the specified sequence or subsequence from the refget store.

        :param identifier: GA4GH sequence identifier (``ga4gh:SQ.…`` or ``SQ.…``)
        :param start: Optional start position (interbase, 0-based)
        :param end: Optional end position (interbase)
        :return: Sequence string
        :raises KeyError: if identifier is not a GA4GH SQ form or not found in store
        :raises DataProxyValidationError: if ``start`` is after ``end``
        """
        # A slice would otherwise return "" for it
        validate_start_le_end(identifier, start, end)
        digest = _extract_digest(identifier)
        if digest is None:
            raise KeyError(identifier)

        self._ensure_sequence_loaded(digest)

        # Full sequence shortcut — no bounds needed
        if start is None and end is None:
            record = self.store.get_sequence(digest)
            if record is None:
                raise KeyError(identifier)
            result = record.decode()
            if result is None:
                raise KeyError(identifier)
            return result

        # Resolve None bounds and clamp to [0, length] for gtars (requires usize).
        # This matches Python slice semantics where negative/out-of-range indices
        # are silently clamped.
        if start is None or end is None or start < 0:
            metadata = self.store.get_sequence_metadata(digest)
            if metadata is None:
                raise KeyError(identifier)
            length: int = metadata.length
            s = max(start or 0, 0)
            e = min(end, length) if end is not None else length
        else:
            s = start
            e = end

        if s >= e:
            return ""

        result = self.store.get_substring(digest, s, e)
        if result is None:
            raise KeyError(identifier)
        return result

    @functools.lru_cache(maxsize=DATAPROXY_CACHE_SIZE)  # noqa: B019
    def get_sequence_length(self, identifier: str) -> int:
        """Return the length of the specified sequence from the refget store.

        :param identifier: GA4GH sequence identifier (``ga4gh:SQ.…`` or ``SQ.…``)
        :return: Sequence length
        :raises KeyError: if identifier is not a GA4GH SQ form or not found in store
        """
        digest = _extract_digest(identifier)
        if digest is None:
            raise KeyError(identifier)
        metadata = self.store.get_sequence_metadata(digest)
        if metadata is None:
            raise KeyError(identifier)
        return metadata.length


class RefgetAliasDataProxy(AliasDataProxy):
    """AliasDataProxy backed by a gtars RefgetStore.

    Resolves any-namespace identifiers to sha512t24u digests via the store's
    sequence alias index, and enumerates all known aliases for a digest via
    the store's reverse-lookup API. The ``sha512t24u:`` and ``ga4gh:SQ.``
    aliases are synthesized from the digest rather than stored on disk.
    """

    def __init__(self, store: RefgetStore) -> None:
        """Initialize RefgetAliasDataProxy.

        :param store: A gtars RefgetStore instance. The store's alias
            sidecars are pulled eagerly for remote stores so subsequent
            lookups do not need extra network / disk I/O for the alias index.
        """
        super().__init__()
        self.store = store
        # pull_aliases is a no-op for local stores but is required for remote
        # stores so namespace alias lookups work from the local cache.
        try:
            store.pull_aliases()  # type: ignore[attr-defined]
        except Exception as exc:
            _logger.debug("pull_aliases failed (ignored): %s", exc)

    @functools.lru_cache(maxsize=DATAPROXY_CACHE_SIZE)  # noqa: B019
    def _resolve_to_digest(self, identifier: str) -> str | None:
        """Resolve an identifier in any supported form to a sha512t24u digest.

        Supported forms:
        - ``ga4gh:SQ.<digest>``, ``SQ.<digest>``, ``sha512t24u:<digest>``
          (handled via :func:`_extract_digest`)
        - ``namespace:alias`` (looked up via ``get_sequence_by_alias``)
        - Bare accession like ``NC_000001.11`` or ``LRG_1t1`` (coerced to
          ``refseq:`` / ``lrg:`` / ... via :func:`coerce_accession_namespace`)

        Returns None if the identifier cannot be resolved.
        """
        digest = _extract_digest(identifier)
        if digest is not None:
            # Digest-prefix branch: confirm the sequence actually exists.
            md = self.store.get_sequence_metadata(digest)
            return digest if md is not None else None

        # Namespaced or bare accession: coerce and look up.
        normalized = identifier
        if ":" not in normalized[1:]:
            try:
                normalized = coerce_accession_namespace(normalized)
            except ValueError:
                return None

        ns, _, alias = normalized.partition(":")
        if not alias:
            return None
        metadata = self.store.get_sequence_metadata_by_alias(ns, alias)
        if metadata is None:
            return None
        return metadata.sha512t24u

    def derive_refget_accession(self, ac: str) -> str | None:
        """Resolve an accession directly without scanning its reverse aliases."""
        if ac is None:
            return None
        digest = self._resolve_to_digest(ac)
        return f"SQ.{digest}" if digest is not None else None

    @functools.lru_cache(maxsize=DATAPROXY_CACHE_SIZE)  # noqa: B019
    def get_aliases(self, identifier: str) -> list[str]:
        """Return all known aliases for the specified sequence.

        Aliases are emitted as ``namespace:alias`` strings. In addition to
        whatever the store's reverse-lookup index returns, two synthesized
        entries are always appended when the sequence exists:

        - ``sha512t24u:<digest>``
        - ``ga4gh:SQ.<digest>``

        :raises KeyError: if the identifier cannot be resolved in the store.
        """
        digest = self._resolve_to_digest(identifier)
        if digest is None:
            raise KeyError(identifier)

        aliases: list[str] = []
        try:
            rev = self.store.get_aliases_for_sequence(digest) or []
        except Exception as exc:
            _logger.debug(
                "get_aliases_for_sequence(%s) raised, treating as empty: %s",
                digest,
                exc,
            )
            rev = []
        for entry in rev:
            if isinstance(entry, tuple) and len(entry) == 2:
                ns, alias = entry
                aliases.append(f"{ns}:{alias}")
            else:
                aliases.append(str(entry))

        aliases.append(f"sha512t24u:{digest}")
        aliases.append(f"ga4gh:SQ.{digest}")
        return aliases


def _create_refget_dataproxy(
    parsed_uri: ParseResult,
    proto: str,
    *,
    disable_healthcheck: bool = False,  # noqa: ARG001
) -> DataProxy:
    """Create a refget-backed DataProxy from a parsed URI.

    Composes a ``RefgetSequenceDataProxy`` and ``RefgetAliasDataProxy`` over
    a single shared ``RefgetStore`` instance, wrapped in a
    ``CompositeDataProxy`` so the result satisfies the full ``DataProxy``
    interface. ``refget+file://`` URIs open a local store directly.
    ``refget+http://`` and ``refget+https://`` URIs open a remote store and
    cache its index, aliases, and fetched sequence payloads locally under
    ``GA4GH_VRS_REFGET_CACHE_DIR`` or the default XDG-style cache root.

    :param parsed_uri: Parsed URI components
    :param proto: Protocol portion of the scheme (e.g. "file", "http", "")
    :param disable_healthcheck: Unused, accepted for interface consistency
    :return: A ``CompositeDataProxy`` wrapping the sequence + alias proxies
    :raises ValueError: if the protocol is not supported
    """
    if proto in ("", "file"):
        from gtars.refget import RefgetStore

        store = RefgetStore.open_local(parsed_uri.path)
    elif proto in ("http", "https"):
        from gtars.refget import RefgetStore

        remote_url = urlunparse(parsed_uri._replace(scheme=proto))
        cache_dir = _resolve_remote_cache_dir(remote_url)
        cache_dir.mkdir(parents=True, exist_ok=True)
        store = RefgetStore.open_remote(str(cache_dir), remote_url)
        # open_remote defers the sequence index, and alias lookups do not load
        # it, so they would find nothing. Metadata lookups load it; the digest
        # is irrelevant.
        store.get_sequence_metadata("")
    else:
        msg = f"Refget URI protocol {proto!r} not supported"
        raise ValueError(msg)

    store.set_quiet(True)
    sequence_proxy = RefgetSequenceDataProxy(store)
    alias_proxy = RefgetAliasDataProxy(store)
    return CompositeDataProxy(sequence_proxy, alias_proxy)
