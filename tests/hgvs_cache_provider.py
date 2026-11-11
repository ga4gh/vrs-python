"""Wrap hgvs provider lookups for offline test recording and replay."""

import functools
from collections.abc import Callable
from typing import Any

import hgvs.dataproviders.uta
import pytest
from hgvs.decorators.lru_cache import _make_key
from hgvs.exceptions import HGVSDataNotAvailableError


def _wrap_cached_method(
    provider: hgvs.dataproviders.uta.UTABase,
    method: Callable[..., Any],
    recording: bool,
) -> Callable[..., Any]:
    """Record and replay "no data" results, and fail on lookups never recorded.

    hgvs caches only returned values. Without this, a lookup that raised
    HGVSDataNotAvailableError would be missing from the recording, and replaying it
    would look the same as a cache miss. A recorded result is stored as the exception
    itself, which hgvs's cache returns like any other value.

    hgvs reports a SeqRepo outage as HGVSDataNotAvailableError too. An outage while
    recording fails tests, so no JSON is written, and any that slips through shows in
    the JSON diff.
    :param provider: Provider whose cache stores lookup results.
    :param method: Cached hgvs lookup method.
    :param recording: Whether live lookups are being recorded.
    :return: Lookup wrapper that preserves recorded missing-data errors.
    """
    name: str = method.__wrapped__.__name__  # hgvs keys cache entries by this name

    @functools.wraps(method)
    def wrapped(*args: object, **kwargs: object) -> Any:  # noqa: ANN401
        """Call the lookup, recording or replaying missing-data errors.

        :param args: Positional lookup arguments.
        :param kwargs: Keyword lookup arguments.
        :return: The cached or live lookup result.
        :raises HGVSDataNotAvailableError: If the lookup found no data.
        """
        try:
            result = method(*args, **kwargs)
        except HGVSDataNotAvailableError as exc:
            if not recording:
                # Replay never queries a data source, so only a cache miss raises.
                # pytest.fail escapes the production handlers that catch this error.
                pytest.fail(
                    f"Missing HGVS recording: {exc}\nRun make record-hgvs-cache."
                )
            # Keep the message: hgvs's Normalizer.is_valid_pos reads it.
            provider.cache[_make_key(name, args, kwargs, False, ())] = (
                HGVSDataNotAvailableError(str(exc))
            )
            provider.cache.sync()
            raise
        if isinstance(result, HGVSDataNotAvailableError):
            # A fresh instance, so tracebacks don't pile up on the cached one.
            raise HGVSDataNotAvailableError(str(result))
        return result

    return wrapped
