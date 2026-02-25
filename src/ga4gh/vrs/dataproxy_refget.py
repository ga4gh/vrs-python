"""Refget/gtars-based implementation of the SequenceDataProxy interface.

Provides a SequenceDataProxy backed by a gtars RefgetStore for high-performance
sequence retrieval by GA4GH digest. Does not provide alias resolution — combine
with an AliasDataProxy via CompositeDataProxy for full DataProxy functionality.

Requires the ``refget`` optional dependency group: ``pip install ga4gh.vrs[refget]``
"""

from __future__ import annotations

import functools
from typing import TYPE_CHECKING
from urllib.parse import ParseResult

from ga4gh.vrs.dataproxy import SequenceDataProxy

if TYPE_CHECKING:
    from gtars.refget import RefgetStore

_GA4GH_SQ_PREFIX = "ga4gh:SQ."
_SQ_PREFIX = "SQ."


def _extract_digest(identifier: str) -> str | None:
    """Strip GA4GH SQ prefix, returning the raw sha512t24u digest.

    Accepts identifiers in ``ga4gh:SQ.<digest>`` or ``SQ.<digest>`` form.

    :param identifier: Sequence identifier, possibly with GA4GH prefix
    :return: Raw digest string, or ``None`` if identifier is not in GA4GH SQ form
    """
    if identifier.startswith(_GA4GH_SQ_PREFIX):
        return identifier[len(_GA4GH_SQ_PREFIX) :]
    if identifier.startswith(_SQ_PREFIX):
        return identifier[len(_SQ_PREFIX) :]
    return None


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

    def get_sequence(
        self, identifier: str, start: int | None = None, end: int | None = None
    ) -> str:
        """Return the specified sequence or subsequence from the refget store.

        :param identifier: GA4GH sequence identifier (``ga4gh:SQ.…`` or ``SQ.…``)
        :param start: Optional start position (interbase, 0-based)
        :param end: Optional end position (interbase)
        :return: Sequence string
        :raises KeyError: if identifier is not a GA4GH SQ form or not found in store
        """
        digest = _extract_digest(identifier)
        if digest is None:
            raise KeyError(identifier)

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
            record = self.store.get_sequence(digest)
            if record is None:
                raise KeyError(identifier)
            length: int = record.metadata.length
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

    @functools.lru_cache(maxsize=512)  # noqa: B019
    def get_sequence_length(self, identifier: str) -> int:
        """Return the length of the specified sequence from the refget store.

        :param identifier: GA4GH sequence identifier (``ga4gh:SQ.…`` or ``SQ.…``)
        :return: Sequence length
        :raises KeyError: if identifier is not a GA4GH SQ form or not found in store
        """
        digest = _extract_digest(identifier)
        if digest is None:
            raise KeyError(identifier)
        record = self.store.get_sequence(digest)
        if record is None:
            raise KeyError(identifier)
        return record.metadata.length


def _create_refget_dataproxy(
    parsed_uri: ParseResult,
    proto: str,
    *,
    disable_healthcheck: bool = False,  # noqa: ARG001
) -> RefgetSequenceDataProxy:
    """Create a RefgetSequenceDataProxy from a parsed URI.

    :param parsed_uri: Parsed URI components
    :param proto: Protocol portion of the scheme (e.g. "file", "")
    :param disable_healthcheck: Unused, accepted for interface consistency
    :return: A RefgetSequenceDataProxy instance
    :raises ValueError: if protocol is not supported
    """
    if proto in ("", "file"):
        from gtars.refget import RefgetStore

        return RefgetSequenceDataProxy(RefgetStore.open_local(parsed_uri.path))
    msg = f"Refget URI protocol {proto!r} not supported"
    raise ValueError(msg)
