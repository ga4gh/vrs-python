"""Fail immediately on live connections or missing HGVS recordings.

pytest.fail raises outside Exception, so production fallback handlers cannot swallow
these failures. VCR replay does not open sockets; psycopg2 needs its own patch because
libpq opens connections in C.
"""

import functools
import socket
from collections.abc import Callable
from typing import Any, NoReturn

import psycopg2
import pytest
from hgvs.exceptions import HGVSDataNotAvailableError


def block_network(mp: pytest.MonkeyPatch) -> None:
    """Block Python socket connections and psycopg2 for offline tests"""

    def forbidden_connection(*_args: object, **_kwargs: object) -> NoReturn:
        pytest.fail("Live network connections are disabled in offline tests.")

    def forbidden_postgres(*_args: object, **_kwargs: object) -> NoReturn:
        pytest.fail("Live PostgreSQL connections are disabled in offline tests.")

    mp.setattr(socket.socket, "connect", forbidden_connection)
    mp.setattr(socket.socket, "connect_ex", forbidden_connection)
    mp.setattr(psycopg2, "connect", forbidden_postgres)


def require_recorded_data(method: Callable[..., Any]) -> Callable[..., Any]:
    """Make a replay-mode provider cache miss fatal, even inside HGVS fallbacks"""

    @functools.wraps(method)
    def wrapped(*args: object, **kwargs: object) -> Any:  # noqa: ANN401
        try:
            return method(*args, **kwargs)
        except HGVSDataNotAvailableError as exc:
            pytest.fail(
                f"Missing HGVS recording: {exc}\nRun make record-hgvs-cache.",
                pytrace=False,
            )

    return wrapped
