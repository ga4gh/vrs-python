import tempfile
import ast
import os
import pytest
import shutil
import sys

from ga4gh.vrs.extras.object_store import Sqlite3MutableMapping


def test_simple(tmp_path):
    db_path = str(tmp_path) + "/test_simple.sqlite3"
    object_store = Sqlite3MutableMapping(db_path)

    kvp = {
        chr(ord("A") + i): i
        for i in range(10)
    }

    assert len(object_store) == 0
    for k, v in kvp.items():
        object_store[k] = v

    assert len(kvp) == len(object_store)
    assert set([k for k, v in kvp.items()]) == set(object_store.keys())

    for k_act, v_act in object_store.items():
        assert kvp[k_act] == v_act

    # Test deletes
    del object_store["A"]
    assert len(object_store) == len(kvp) - 1
    with pytest.raises(KeyError):
        del object_store["A"]
    while len(object_store) > 0:
        del object_store[list(object_store.keys())[0]]
    assert len(object_store) == 0


def test_complex(tmp_path):
    db_path = str(tmp_path) + "/test_complex.sqlite3"

    kvp = {
        "A": "A-value",
        "B": {
            "B-1": "B-1-value",
            "B-2": {
                "B-2-1": [
                    "B-2-1-1",
                    "B-2-1-2",
                    "B-2-1-3"
                ],
                "B-3": 12345
            }
        }
    }

    object_store = Sqlite3MutableMapping(db_path)
    assert len(object_store) == 0

    for k, v in kvp.items():
        object_store[k] = v

    assert len(kvp) == len(object_store)
    assert set([k for k, v in kvp.items()]) == set(object_store.keys())

    for k_act, v_act in object_store.items():
        assert kvp[k_act] == v_act

    object_store.close()


def test_classes(tmp_path):
    db_path = str(tmp_path) + "/test_complex.sqlite3"

    class TestClass(object):
        def __init__(self, id):
            self.id = id
            self.A = "A"
            self.B = ["B1", "B2", 3]
            self.C = {"C1": "C1-value"}

        def somefunction(self, arg):
            return f"{self.id}-{arg}"

    object_store = Sqlite3MutableMapping(db_path)

    val1 = TestClass("val1")
    val2 = TestClass("val2")
    object_store["val1"] = val1
    object_store["val2"] = val2

    # Test retrieval of custom object contents
    assert len(object_store) == 2
    assert object_store["val1"].id == "val1"
    assert object_store["val2"].id == "val2"
    assert object_store["val1"].somefunction("X") == "val1-X"

    # Test deletes of custom objects
    del object_store["val1"]
    assert len(object_store) == 1
    with pytest.raises(KeyError):
        del object_store["val1"]

    del object_store["val2"]
    assert len(object_store) == 0

    object_store.close()


# This version verifies that setting autocommit=False may
# lose data is .close or .commit is not explicitly called
# def test_no_commit():
#     tmpdir = tempfile.mkdtemp()
#     db_path = tmpdir + "/test_commit.sqlite3"
#     object_store = Sqlite3MutableMapping(db_path, autocommit=False)
#     value_count = int(1e5)
#     for i in range(value_count):
#         object_store[f"key{i}"] = f"value{i}"
#     object_store.db.close()
#     # See if the stuff is still there
#     object_store = Sqlite3MutableMapping(db_path)
#     for i in range(value_count):
#         assert object_store[f"key{i}"] == f"value{i}"

def test_commit(tmp_path):
    db_path = str(tmp_path) + "/test_commit.sqlite3"
    object_store = Sqlite3MutableMapping(db_path, autocommit=False)

    value_count = int(1e4)
    for i in range(value_count):
        object_store[f"key{i}"] = f"value{i}"

    object_store.close()

    # See if the stuff is still there
    object_store = Sqlite3MutableMapping(db_path)
    for i in range(value_count):
        assert object_store[f"key{i}"] == f"value{i}"


def test_contextmanager(tmp_path):
    db_path = str(tmp_path) + "/test_contextmanager.sqlite3"
    value_count = int(1e4)
    with Sqlite3MutableMapping(db_path, autocommit=False) as object_store:
        for i in range(value_count):
            object_store[f"key{i}"] = f"value{i}"

    with Sqlite3MutableMapping(db_path, autocommit=False) as object_store:
        # See if the stuff is still there
        for i in range(value_count):
            assert object_store[f"key{i}"] == f"value{i}"
