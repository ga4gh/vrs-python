"""Refget data proxies over a real gtars RefgetStore, opened locally and over HTTP.

The store is built from the sequences and aliases of the SeqRepo in
``tests/data/seqrepo/latest``, so both backends serve the same data.
"""

import shutil
import sqlite3
import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path

import pytest

from ga4gh.vrs.dataproxy import DataProxy, create_dataproxy

SEQREPO_DIR = Path(__file__).parent / "data" / "seqrepo" / "latest"
# SeqRepo alias namespace -> RefgetStore alias namespace
NAMESPACES = {"NCBI": "refseq", "GRCh38": "GRCh38"}
VHL_DIGEST = "v_QTc1p-MUYdgrRv4LMT6ByXIOsdw3C_"  # NM_000551.3


@pytest.fixture(scope="session")
def refget_store_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    refget = pytest.importorskip("gtars.refget")
    work_dir = tmp_path_factory.mktemp("refget")
    store_dir = work_dir / "store"
    store = refget.RefgetStore.on_disk(str(store_dir))
    store.set_quiet(True)
    # SeqRepo names each FASTA record by its sha512t24u digest. Import copies:
    # gtars writes a digest index next to each FASTA it reads.
    for fasta in sorted((SEQREPO_DIR / "sequences").rglob("*.fa.bgz")):
        copy = work_dir / fasta.name
        shutil.copyfile(fasta, copy)
        store.add_sequence_collection_from_fasta(str(copy))
    digests = {metadata.sha512t24u for metadata in store.list_sequences()}
    with sqlite3.connect(SEQREPO_DIR / "aliases.sqlite3") as aliases:
        for seq_id, namespace, alias in aliases.execute(
            "SELECT seq_id, namespace, alias FROM seqalias WHERE is_current = 1"
        ):
            if namespace in NAMESPACES and seq_id in digests:
                store.add_sequence_alias(NAMESPACES[namespace], alias, seq_id)
    store.write()
    return store_dir


@pytest.fixture(scope="session")
def refget_store_url(refget_store_dir: Path) -> Iterator[str]:
    # A separate process: gtars holds the GIL during its HTTP requests, so an
    # in-process server thread could never answer them
    server = subprocess.Popen(  # noqa: S603
        [
            sys.executable,
            "-u",
            "-m",
            "http.server",
            "0",
            "--bind",
            "127.0.0.1",
            "--directory",
            str(refget_store_dir),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    try:
        # "Serving HTTP on 127.0.0.1 port 49482 (http://127.0.0.1:49482/) ..."
        port = int(server.stdout.readline().split(" port ")[1].split()[0])
        yield f"http://127.0.0.1:{port}"
    finally:
        server.terminate()
        server.wait()


@pytest.fixture(params=["local", "remote"])
def refget_dataproxy(
    request: pytest.FixtureRequest,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    refget_store_dir: Path,
) -> DataProxy:
    if request.param == "local":
        return create_dataproxy(f"refget+file://{refget_store_dir}")
    monkeypatch.setenv("GA4GH_VRS_REFGET_CACHE_DIR", str(tmp_path / "cache"))
    return create_dataproxy(f"refget+{request.getfixturevalue('refget_store_url')}")


@pytest.mark.parametrize(
    "identifier", ["NM_000551.3", "refseq:NM_000551.3", f"ga4gh:SQ.{VHL_DIGEST}"]
)
def test_resolves_identifier(refget_dataproxy: DataProxy, identifier: str) -> None:
    assert refget_dataproxy.derive_refget_accession(identifier) == f"SQ.{VHL_DIGEST}"


def test_matches_seqrepo(refget_dataproxy: DataProxy) -> None:
    seqrepo_dataproxy = create_dataproxy(f"seqrepo+file://{SEQREPO_DIR}")
    for identifier in ("NM_000551.3", "NC_000019.10", "GRCh38:chr19"):
        assert refget_dataproxy.get_sequence_length(identifier) == (
            seqrepo_dataproxy.get_sequence_length(identifier)
        )
        assert refget_dataproxy.get_sequence(identifier, 1000, 1100) == (
            seqrepo_dataproxy.get_sequence(identifier, 1000, 1100)
        )
    assert "refseq:NM_000551.3" in refget_dataproxy.get_aliases("NM_000551.3")
