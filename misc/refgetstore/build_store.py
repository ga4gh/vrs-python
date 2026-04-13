#!/usr/bin/env python
"""Build a custom gtars RefgetStore with vrs-python alias parity.

Ingests NCBI genomic FASTA files into a local RefgetStore and augments it with
per-sequence and per-collection aliases parsed from the matching NCBI assembly
reports. See misc/refgetstore/README.md for details.
"""

from __future__ import annotations

import argparse
import logging
import sys
import tomllib
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

from gtars.refget import RefgetStore

logger = logging.getLogger("build_store")

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = Path(__file__).parent / "assemblies.toml"
DEFAULT_STORE = Path(__file__).parent / "store"
DEFAULT_DOWNLOADS = Path(__file__).parent / "downloads"

SQ_PREFIX = "SQ."
NA_VALUES = {"na", "", "<NA>"}


@dataclass
class AssemblyConfig:
    namespace: str
    fasta_url: str
    report_url: str
    load_fasta: bool = True
    fasta_path: str | None = None


@dataclass
class AssemblyReport:
    """Parsed NCBI assembly_report.txt."""

    refseq_assembly_accession: str | None
    genbank_assembly_accession: str | None
    rows: list[dict[str, str]] = field(default_factory=list)


@dataclass
class AssemblyStats:
    namespace: str
    collection_digest: str | None = None
    rows_total: int = 0
    rows_resolved: int = 0
    rows_skipped_non_refseq: int = 0
    sequence_aliases_added: int = 0
    sequence_aliases_skipped: int = 0
    collection_aliases_added: int = 0
    warnings: int = 0


def load_config(path: Path) -> list[AssemblyConfig]:
    with path.open("rb") as fh:
        data = tomllib.load(fh)
    return [AssemblyConfig(**entry) for entry in data.get("assembly", [])]


def ensure_download(url: str, target: Path, force: bool) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not force:
        logger.info("using cached %s", target)
        return target
    logger.info("downloading %s -> %s", url, target)
    tmp = target.with_suffix(target.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)  # noqa: S310
    tmp.replace(target)
    return target


def resolve_fasta_source(
    entry: AssemblyConfig, download_dir: Path, force: bool
) -> Path | None:
    if not entry.load_fasta:
        return None
    if entry.fasta_path:
        local = (REPO_ROOT / entry.fasta_path).resolve()
        if local.exists():
            logger.info("using local fasta override %s", local)
            return local
        logger.warning(
            "fasta_path %s not found; falling back to download", local
        )
    return ensure_download(
        entry.fasta_url,
        download_dir / Path(entry.fasta_url).name,
        force,
    )


def resolve_report_source(
    entry: AssemblyConfig, download_dir: Path, force: bool
) -> Path:
    return ensure_download(
        entry.report_url,
        download_dir / f"{entry.namespace}.assembly_report.txt",
        force,
    )


def parse_assembly_report(path: Path) -> AssemblyReport:
    """Parse an NCBI assembly_report.txt into a structured object.

    Header lines begin with ``#`` and include ``GenBank assembly accession`` /
    ``RefSeq assembly accession`` entries. The final ``#``-prefixed line is the
    column header for the data rows (a tab-separated table).
    """
    refseq_accn: str | None = None
    genbank_accn: str | None = None
    columns: list[str] | None = None
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            if line.startswith("#"):
                stripped = line.lstrip("#").strip()
                if stripped.startswith("GenBank assembly accession:"):
                    genbank_accn = stripped.split(":", 1)[1].strip() or None
                elif stripped.startswith("RefSeq assembly accession:"):
                    refseq_accn = stripped.split(":", 1)[1].strip() or None
                elif "\t" in stripped and "Sequence-Name" in stripped:
                    columns = stripped.split("\t")
                continue
            if columns is None:
                continue
            parts = line.split("\t")
            if len(parts) != len(columns):
                continue
            rows.append(dict(zip(columns, parts)))
    return AssemblyReport(
        refseq_assembly_accession=refseq_accn,
        genbank_assembly_accession=genbank_accn,
        rows=rows,
    )


def build_name_to_digest_map(
    store: RefgetStore, collection_digest: str
) -> dict[str, str]:
    store.load_collection(collection_digest)
    level2 = store.get_collection_level2(collection_digest)
    names = level2["names"]
    sequences = level2["sequences"]
    out: dict[str, str] = {}
    for name, seq_id in zip(names, sequences):
        digest = seq_id[len(SQ_PREFIX) :] if seq_id.startswith(SQ_PREFIX) else seq_id
        out[name] = digest
    return out


def build_global_name_to_digest_map(store: RefgetStore) -> dict[str, str]:
    """Fallback map built by scanning every sequence in the store.

    Used when `load_fasta=false` so we can still resolve refseq accessions to
    a digest without having just ingested this assembly's collection. Names
    are globally unique-ish (seqrepo jungle store: 1224 collisions in 580k);
    we keep the first digest we see for each name.
    """
    out: dict[str, str] = {}
    for record in store.iter_sequences():
        md = record.metadata
        out.setdefault(md.name, md.sha512t24u)
    return out


def add_unique_alias(
    store: RefgetStore,
    namespace: str,
    alias: str,
    digest: str,
    stats: AssemblyStats,
) -> None:
    existing = store.get_sequence_by_alias(namespace, alias)
    if existing is not None:
        stats.sequence_aliases_skipped += 1
        return
    store.add_sequence_alias(namespace, alias, digest)
    stats.sequence_aliases_added += 1


def add_aliases_for_row(
    store: RefgetStore,
    row: dict[str, str],
    namespace: str,
    name_to_digest: dict[str, str],
    stats: AssemblyStats,
) -> None:
    refseq_ac = row.get("RefSeq-Accn", "").strip()
    ucsc_name = row.get("UCSC-style-name", "").strip()
    genbank_ac = row.get("GenBank-Accn", "").strip()

    # Scope: we only alias sequences that are actually in the RefSeq dataset.
    # Rows with RefSeq-Accn=na are GenBank-only contigs that the NCBI
    # _genomic.fna.gz (RefSeq view) deliberately omits, so there is no
    # sequence data to alias to. Skip cleanly without a warning.
    if refseq_ac.lower() in NA_VALUES:
        stats.rows_skipped_non_refseq += 1
        return

    digest = name_to_digest.get(refseq_ac)
    if digest is None:
        stats.warnings += 1
        logger.warning(
            "no digest for %s in namespace %s (FASTA does not contain it)",
            refseq_ac,
            namespace,
        )
        return

    stats.rows_resolved += 1
    if ucsc_name and ucsc_name.lower() not in NA_VALUES:
        add_unique_alias(store, namespace, ucsc_name, digest, stats)
    if genbank_ac and genbank_ac.lower() not in NA_VALUES:
        add_unique_alias(store, namespace, genbank_ac, digest, stats)
        add_unique_alias(store, "insdc", genbank_ac, digest, stats)
    add_unique_alias(store, "refseq", refseq_ac, digest, stats)


def add_collection_aliases(
    store: RefgetStore,
    report: AssemblyReport,
    collection_digest: str,
    stats: AssemblyStats,
) -> None:
    if report.refseq_assembly_accession:
        store.add_collection_alias(
            "refseq", report.refseq_assembly_accession, collection_digest
        )
        stats.collection_aliases_added += 1
    if report.genbank_assembly_accession:
        store.add_collection_alias(
            "insdc", report.genbank_assembly_accession, collection_digest
        )
        stats.collection_aliases_added += 1


def process_assembly(
    store: RefgetStore,
    entry: AssemblyConfig,
    download_dir: Path,
    force_download: bool,
    global_name_map_cache: dict[str, dict[str, str]],
) -> AssemblyStats:
    stats = AssemblyStats(namespace=entry.namespace)
    logger.info("=== %s ===", entry.namespace)
    fasta_path = resolve_fasta_source(entry, download_dir, force_download)
    report_path = resolve_report_source(entry, download_dir, force_download)

    name_to_digest: dict[str, str]
    coll_digest: str | None = None
    if fasta_path is not None:
        logger.info("ingesting %s", fasta_path)
        coll_meta, was_new = store.add_sequence_collection_from_fasta(
            str(fasta_path)
        )
        coll_digest = coll_meta.digest
        stats.collection_digest = coll_digest
        logger.info(
            "collection %s (%s, %d sequences)",
            coll_digest,
            "new" if was_new else "existing",
            coll_meta.n_sequences,
        )
        name_to_digest = build_name_to_digest_map(store, coll_digest)
    else:
        logger.info("load_fasta=false; using global name map")
        cached = global_name_map_cache.get("_global")
        if cached is None:
            cached = build_global_name_to_digest_map(store)
            global_name_map_cache["_global"] = cached
        name_to_digest = cached

    logger.info("parsing %s", report_path)
    report = parse_assembly_report(report_path)
    stats.rows_total = len(report.rows)
    logger.info(
        "report: %d rows, refseq=%s genbank=%s",
        stats.rows_total,
        report.refseq_assembly_accession,
        report.genbank_assembly_accession,
    )

    for row in report.rows:
        add_aliases_for_row(store, row, entry.namespace, name_to_digest, stats)

    if coll_digest is not None:
        add_collection_aliases(store, report, coll_digest, stats)

    return stats


def print_summary(all_stats: list[AssemblyStats], store: RefgetStore) -> None:
    print()
    print("Summary")
    print("-------")
    for s in all_stats:
        print(
            f"  {s.namespace:14s}  rows={s.rows_total:4d} "
            f"resolved={s.rows_resolved:4d} "
            f"non_refseq={s.rows_skipped_non_refseq:3d} "
            f"seq_aliases_added={s.sequence_aliases_added:4d} "
            f"skipped={s.sequence_aliases_skipped:4d} "
            f"coll_aliases_added={s.collection_aliases_added:2d} "
            f"warnings={s.warnings:3d}"
        )
    print()
    print("Store stats:", store.stats())
    print(
        "Sequence alias namespaces:",
        sorted(store.list_sequence_alias_namespaces()),
    )
    print(
        "Collection alias namespaces:",
        sorted(store.list_collection_alias_namespaces()),
    )


def iter_selected(
    entries: list[AssemblyConfig], selected: str | None
) -> Iterator[AssemblyConfig]:
    if selected is None:
        yield from entries
        return
    found = False
    for entry in entries:
        if entry.namespace == selected:
            found = True
            yield entry
    if not found:
        raise SystemExit(
            f"no assembly with namespace {selected!r} in config"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--store-dir", type=Path, default=DEFAULT_STORE)
    parser.add_argument("--download-dir", type=Path, default=DEFAULT_DOWNLOADS)
    parser.add_argument("--assembly", type=str, default=None)
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    entries = load_config(args.config)
    logger.info("loaded %d assembly entries from %s", len(entries), args.config)

    args.store_dir.mkdir(parents=True, exist_ok=True)
    args.download_dir.mkdir(parents=True, exist_ok=True)

    store = RefgetStore.on_disk(str(args.store_dir))
    logger.info("opened store at %s (mode=%s)", args.store_dir, store.storage_mode)

    all_stats: list[AssemblyStats] = []
    global_name_map_cache: dict[str, dict[str, str]] = {}
    for entry in iter_selected(entries, args.assembly):
        stats = process_assembly(
            store,
            entry,
            args.download_dir,
            args.force_download,
            global_name_map_cache,
        )
        all_stats.append(stats)

    logger.info("persisting store to %s", args.store_dir)
    store.write()
    print_summary(all_stats, store)


if __name__ == "__main__":
    sys.exit(main())
