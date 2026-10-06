"""Convert between hgvs's pickled data provider cache and the JSON file tests commit.

hgvs reads and writes its cache (`hgvs.dataproviders.uta.connect(mode=..., cache=...)`)
as a pickle, which can run code when loaded and can't be reviewed in a diff. Tests
instead commit `tests/data/hgvs_cache.json` and generate the pickle from it at startup.

Each cache entry is stored as ``{"func": ..., "args": [...], "value": ...}``. Values are
strings, numbers, None, lists, and UTA rows (`psycopg2.extras.DictRow`), which are
stored as ``{"psycopg2_dictrow": {column: value, ...}}`` in column order.

Usage, to convert a cache hgvs just recorded:
``python tests/hgvs_cache_json.py <recorded.pkl> <hgvs_cache.json>``
"""

import json
import pickle
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Any

from hgvs.decorators.lru_cache import _make_key
from psycopg2.extras import DictRow

# hgvs builds each cache key as (*args, "__func__", function name)
_FUNC_MARKER = "__func__"
_DICTROW_TAG = "psycopg2_dictrow"
_JSON_SCALARS = (str, int, float, bool, type(None))


def _encode_value(value: Any) -> Any:  # noqa: ANN401
    """Convert a cached hgvs result to a JSON value"""
    if isinstance(value, DictRow):
        return {_DICTROW_TAG: dict(value.items())}
    if isinstance(value, list):
        return [_encode_value(v) for v in value]
    if isinstance(value, _JSON_SCALARS):
        return value
    msg = f"can't store a {type(value).__name__} in the hgvs JSON cache: {value!r}"
    raise TypeError(msg)


def _decode_value(value: Any) -> Any:  # noqa: ANN401
    """Convert a JSON value back to the cached hgvs result, rebuilding UTA rows"""
    if isinstance(value, dict):
        if value.keys() != {_DICTROW_TAG}:
            msg = f'expected {{"{_DICTROW_TAG}": {{...}}}}, got {value!r}'
            raise ValueError(msg)
        columns: dict[str, Any] = value[_DICTROW_TAG]
        row = DictRow.__new__(DictRow)
        index = OrderedDict((column, i) for i, column in enumerate(columns))
        row.__setstate__((list(columns.values()), index))
        return row
    if isinstance(value, list):
        return [_decode_value(v) for v in value]
    return value


def pickle_to_json(pickle_file: Path, json_file: Path) -> None:
    """Write the hgvs cache pickle `pickle_file` as JSON to `json_file`

    Only use this on a cache hgvs just recorded: loading a pickle can run code.
    """
    with pickle_file.open("rb") as f:
        cache: dict = pickle.load(f)  # noqa: S301  # hgvs cache key -> cached result
    entries = []
    for key in cache:
        *args, marker, func = key
        if marker != _FUNC_MARKER:
            msg = f"unexpected hgvs cache key layout: {key!r}"
            raise ValueError(msg)
        entries.append({"func": func, "args": args, "value": _encode_value(cache[key])})
    entries.sort(key=lambda e: (e["func"], json.dumps(e["args"])))
    json_file.write_text(json.dumps(entries, indent=2) + "\n")


def json_to_pickle(json_file: Path, pickle_file: Path) -> None:
    """Write the JSON hgvs cache `json_file` as the pickle hgvs loads to `pickle_file`"""
    entries: list[dict[str, Any]] = json.loads(json_file.read_text())
    cache = {}
    for e in entries:
        try:
            value = _decode_value(e["value"])
        except ValueError as err:
            msg = f"{json_file}: {e['func']} {e['args']}: {err}"
            raise ValueError(msg) from err
        cache[_make_key(e["func"], tuple(e["args"]), {}, False, ())] = value
    with pickle_file.open("wb") as f:
        pickle.dump(cache, f)


if __name__ == "__main__":
    pickle_to_json(Path(sys.argv[1]), Path(sys.argv[2]))
