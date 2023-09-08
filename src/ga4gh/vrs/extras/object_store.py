from collections.abc import MutableMapping
from typing import Any, Union
from threading import Lock

import sqlite3
import dill


class Sqlite3MutableMapping(MutableMapping):
    """
    Class that can be used like a Python dictionary but that uses a sqlite3 database
    as the storage. Can also be opened as a contextmanager.

    If not used as a contextmanager, user must call commit and/or close.
    """

    def __init__(
            self,
            sqlite3_db: Union[str, sqlite3.Connection],
            autocommit=True
    ):
        """
        Connect to the sqlite3 database specified by an existing sqlite3.Connection
        or a connection string.

        - autocommit: if False, disables commit after every setitem/delitem.
                Significant performance implication (>10X speedup)
        """
        if isinstance(sqlite3_db, str):
            sqlite3_db = sqlite3.connect(
                sqlite3_db,
                check_same_thread=True)
        self.db = sqlite3_db
        self.autocommit = autocommit
        self._closed_lock = Lock()
        self._closed = False
        self._create_schema()

    def _create_schema(self):
        cur = self.db.cursor()
        try:
            cur.execute(
                "create table if not exists mapping "
                "(key text, value blob)")
            cur.execute(
                "create unique index if not exists mapping_key_idx "
                "on mapping (key)")
            self.commit()
        finally:
            cur.close()

    def __del__(self):
        self.close()

    def __delitem__(self, key: Any) -> None:
        # Raise KeyError
        self[key]
        # Delete if found
        cur = self.db.cursor()
        try:
            cur.execute(
                "delete from mapping where key = ?",
                (key,))
            if self.autocommit:
                self.commit()
        finally:
            cur.close()

    def __setitem__(self, key: Any, value: Any) -> None:
        cur = self.db.cursor()
        try:
            ser = dill.dumps(value)
            cur.execute(
                "insert or replace into mapping(key, value) "
                "values (?, ?)",
                (key, sqlite3.Binary(ser)))
            if self.autocommit:
                self.commit()
        finally:
            cur.close()

    def __getitem__(self, key: Any) -> Any:
        cur = self.db.cursor()
        try:
            rows = cur.execute(
                "select value from mapping where key = ?",
                (key,))
            row0 = next(rows)
            if row0:
                des = dill.loads(row0[0])
                return des
        except StopIteration:
            raise KeyError("Key not found: " + str(key))
        finally:
            cur.close()

    def __iter__(self):
        cur = self.db.cursor()
        try:
            rows = cur.execute("select key from mapping")
            for row in rows:
                yield row[0]
        finally:
            cur.close()

    def __len__(self):
        cur = self.db.cursor()
        try:
            rows = cur.execute("select count(*) from mapping")
            ct = list(rows)[0][0]
            return ct
        finally:
            cur.close()

    def commit(self):
        self.db.commit()

    def close(self):
        with self._closed_lock:
            if not self._closed:
                self.commit()
                self.db.close()
            self._closed = True

    def __enter__(self):
        self.db.__enter__()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.db.__exit__(exc_type, exc_value, traceback)
