"""Test recording hgvs lookups, which CI otherwise never runs (see conftest.py)"""

from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import Any, NoReturn

import pytest
from conftest import _wrap_cached_method
from hgvs.decorators.lru_cache import LEARN, RUN, lru_cache
from hgvs.exceptions import HGVSDataNotAvailableError
from hgvs.utils.PersistentDict import PersistentDict
from hgvs_cache_json import json_to_pickle, pickle_to_json


def test_not_found_survives_record_and_replay(tmp_path: Path) -> None:
    def get_tx_identity_info(ac: str) -> NoReturn:
        msg = f"No data for {ac}"
        raise HGVSDataNotAvailableError(msg)

    def wrap(cache: PersistentDict, mode: int) -> Callable[[str], Any]:
        cached = lru_cache(maxsize=None, mode=mode, cache=cache)(get_tx_identity_info)
        return _wrap_cached_method(SimpleNamespace(cache=cache), cached, mode == LEARN)

    with pytest.raises(HGVSDataNotAvailableError):
        wrap(PersistentDict(str(tmp_path / "rec.pkl")), LEARN)("NM_X")
    pickle_to_json(tmp_path / "rec.pkl", tmp_path / "c.json")
    json_to_pickle(tmp_path / "c.json", tmp_path / "rep.pkl")
    # An unrecorded lookup would fail the test here rather than raise this error.
    with pytest.raises(HGVSDataNotAvailableError, match="^No data for NM_X$"):
        wrap(PersistentDict(str(tmp_path / "rep.pkl"), flag="r"), RUN)("NM_X")
