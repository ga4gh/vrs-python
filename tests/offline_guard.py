"""Keep tests from going past their recorded data.

Tests replay recorded data: VCR cassettes for HTTP requests, and the hgvs cache for
hgvs's UTA queries and sequence fetches. `conftest.py` installs an `OfflineGuard` for
the session that:

- with `--vcr-record=none` (as in CI), blocks network connections, both through Python
  sockets (all HTTP clients) and `psycopg2.connect` (libpq connects in C). VCR serves
  recorded requests without opening a connection, so any connection is a request
  missing from a cassette, or from a test not marked `@pytest.mark.vcr`.
- with the hgvs cache in run mode (the default), records cache misses, which
  `HgvsTools` would otherwise catch and only log.

Each violation is recorded as well as raised, because library code catches some
errors, so a raised error alone could be swallowed. `conftest.py` fails the test during
which a violation was recorded.
"""

import functools
import socket
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any, NoReturn

import hgvs.dataproviders.uta
import psycopg2
import pytest
from hgvs.exceptions import HGVSDataNotAvailableError

# Socket families that reach the network; others, like Unix sockets, are local IPC
_NETWORK_FAMILIES = (socket.AF_INET, socket.AF_INET6)


class OfflineError(AssertionError):
    """A test opened a network connection while tests replay recorded data"""


def _describe_call(name: str, args: tuple, kwargs: dict[str, Any]) -> str:
    params = [repr(a) for a in args] + [f"{k}={v!r}" for k, v in kwargs.items()]
    return f"{name}({', '.join(params)})"


@dataclass
class OfflineGuard:
    """Blocks network connections and records hgvs cache misses while installed"""

    violations: list[str] = field(default_factory=list)
    blocking_network: bool = False
    recording_hgvs_cache_misses: bool = False

    def block_network(self, mp: pytest.MonkeyPatch) -> None:
        """Make network connections raise `OfflineError`"""
        for name in ("connect", "connect_ex"):
            mp.setattr(socket.socket, name, self._blocked_socket_method(name))
        mp.setattr(psycopg2, "connect", self._blocked("psycopg2.connect"))
        self.blocking_network = True

    def record_hgvs_cache_misses(
        self, mp: pytest.MonkeyPatch, provider: hgvs.dataproviders.uta.UTABase
    ) -> None:
        """Record cache misses of `provider`, which must use the hgvs cache in run mode

        hgvs wraps each of its lookups in its own `lru_cache`, stored as an instance
        attribute. In run mode that wrapper raises `HGVSDataNotAvailableError` on a miss
        without calling the lookup, so any such error coming out of these methods is a
        cache miss.
        """
        cached_lookups = {
            name: method
            for name, method in vars(provider).items()
            if hasattr(method, "cache_info")
        }
        if not cached_lookups:
            msg = "found no cached hgvs data provider lookups to guard"
            raise RuntimeError(msg)
        for name, method in cached_lookups.items():
            mp.setattr(provider, name, self._recording_misses(name, method))
        self.recording_hgvs_cache_misses = True

    def _block(self, call: str) -> NoReturn:
        self.violations.append(f"blocked: {call}")
        msg = f"{call} opens a network connection while tests replay recorded data"
        raise OfflineError(msg)

    def _blocked(self, name: str) -> Callable[..., NoReturn]:
        def blocked(*args: object, **kwargs: object) -> NoReturn:
            self._block(_describe_call(name, args, kwargs))

        return blocked

    def _blocked_socket_method(self, name: str) -> Callable[..., Any]:
        original = getattr(socket.socket, name)

        def blocked_socket_method(sock: socket.socket, *args: object) -> Any:  # noqa: ANN401
            if sock.family in _NETWORK_FAMILIES:
                self._block(_describe_call(f"socket.{name}", args, {}))
            return original(sock, *args)

        return blocked_socket_method

    def _recording_misses(
        self, name: str, method: Callable[..., Any]
    ) -> Callable[..., Any]:
        @functools.wraps(method)
        def recording_misses(*args: object, **kwargs: object) -> Any:  # noqa: ANN401
            try:
                return method(*args, **kwargs)
            except HGVSDataNotAvailableError:
                self.violations.append(
                    f"hgvs cache miss: {_describe_call(name, args, kwargs)}"
                )
                raise

        return recording_misses
