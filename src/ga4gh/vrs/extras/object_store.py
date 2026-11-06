"""Store pickled Python objects in a SQLite-backed mutable mapping."""

import pickle
import sqlite3
from collections.abc import Iterator, MutableMapping
from threading import Lock
from types import TracebackType
from typing import Self


class Sqlite3MutableMapping(MutableMapping[str, object]):
    """Store a mutable mapping in SQLite.

    The context manager commits or rolls back pending transactions. With
    autocommit enabled, earlier writes have already been committed and cannot
    be rolled back by the context manager. Call ``close`` to
    release the connection. With autocommit disabled, call ``commit`` to persist
    changes outside a context manager.
    """

    def __init__(
        self, sqlite3_db: str | sqlite3.Connection, autocommit: bool = True
    ) -> None:
        """Open a database connection and initialize the mapping schema.

        :param sqlite3_db: Database path or existing connection. Schema initialization
            commits pending work on an existing connection; close also closes it.
        :param autocommit: Commit after each insertion or deletion when true.
        """
        if isinstance(sqlite3_db, str):
            sqlite3_db = sqlite3.connect(sqlite3_db, check_same_thread=True)
        self.db = sqlite3_db
        self.autocommit = autocommit
        self._closed_lock = Lock()
        self._closed = False
        self._create_schema()

    def _create_schema(self) -> None:
        """Create the mapping table and unique key index."""
        cur = self.db.cursor()
        try:
            cur.execute("create table if not exists mapping (key text, value blob)")
            cur.execute(
                "create unique index if not exists mapping_key_idx on mapping (key)"
            )
            self.commit()
        finally:
            cur.close()

    def __del__(self) -> None:
        """Commit pending changes and close the connection."""
        self.close()

    def __delitem__(self, key: str) -> None:
        """Delete a stored value.

        :param key: Key to remove.
        :raises KeyError: If the key is absent.
        """
        self[key]
        cur = self.db.cursor()
        try:
            cur.execute("delete from mapping where key = ?", (key,))
            if self.autocommit:
                self.commit()
        finally:
            cur.close()

    def __setitem__(self, key: str, value: object) -> None:
        """Pickle and store a value, replacing any existing entry.

        :param key: Storage key.
        :param value: Picklable value to store.
        """
        cur = self.db.cursor()
        try:
            ser = pickle.dumps(value)
            cur.execute(
                "insert or replace into mapping(key, value) values (?, ?)",
                (key, sqlite3.Binary(ser)),
            )
            if self.autocommit:
                self.commit()
        finally:
            cur.close()

    def __getitem__(self, key: str) -> object:
        """Load and unpickle a stored value.

        :param key: Storage key.
        :return: Stored Python object.
        :raises KeyError: If the key is absent.
        """
        cur = self.db.cursor()
        try:
            rows = cur.execute("select value from mapping where key = ?", (key,))
            row0 = next(rows)
            if row0:
                return pickle.loads(row0[0])  # noqa: S301
        except StopIteration as e:
            raise KeyError("Key not found: " + str(key)) from e
        finally:
            cur.close()

    def __iter__(self) -> Iterator[str]:
        """Iterate over stored keys.

        :return: Iterator yielding each key.
        """
        cur = self.db.cursor()
        try:
            rows = cur.execute("select key from mapping")
            for row in rows:
                yield row[0]
        finally:
            cur.close()

    def __len__(self) -> int:
        """Count stored entries.

        :return: Number of entries.
        """
        cur = self.db.cursor()
        try:
            rows = cur.execute("select count(*) from mapping")
            return list(rows)[0][0]  # noqa: RUF015
        finally:
            cur.close()

    def commit(self) -> None:
        """Commit pending changes to the database."""
        self.db.commit()

    def close(self) -> None:
        """Commit pending changes and close the connection once."""
        with self._closed_lock:
            if not self._closed:
                self.commit()
                self.db.close()
            self._closed = True

    def __enter__(self) -> Self:
        """Enter the database transaction context.

        :return: This mapping.
        """
        self.db.__enter__()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Commit a successful transaction or roll back a failed transaction.

        :param exc_type: Exception class, if the context raised.
        :param exc_value: Exception instance, if the context raised.
        :param traceback: Traceback for the exception, if any.
        """
        self.db.__exit__(exc_type, exc_value, traceback)
