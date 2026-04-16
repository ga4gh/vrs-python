#!/usr/bin/env python
"""Diagnose Ensembl digest mismatches between seqrepo and refgetstore.

A cross-check of the ``ensembl`` namespace found 116 aliases where seqrepo's
stored ``seq_id`` differs from refgetstore's ``sha512t24u`` digest. This
script investigates the root cause by:

    Phase A — Enumerate all mismatches between the two alias tables.
    Phase B — For each mismatch, fetch the full sequence from both stores
              and compare byte-by-byte.
    Phase C — Classify the root cause (``*`` normalization, casing, etc.)
              and determine which digest is correct per GA4GH VRS spec.
    Phase D — Summarize findings and write a JSON report.

Usage:

    SEQREPO_ROOT_DIR=/path/to/seqrepo/2024-12-20 \\
        uv run python misc/refgetstore/diagnose_ensembl_mismatches.py \\
            [--refget-store misc/refgetstore/store] \\
            [--report misc/refgetstore/diagnosis_report.json]

Exits 0 on success (report written), 1 on setup error.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

from biocommons.seqrepo import SeqRepo
from bioutils.digests import seq_seqhash
from ga4gh.core.digests import sha512t24u
from gtars.refget import RefgetStore

HERE = Path(__file__).resolve().parent
DEFAULT_STORE_PATH = HERE / "store"


@dataclass
class MismatchRecord:
    """A single ensembl alias where seqrepo and refgetstore digests differ."""

    accession: str
    seqrepo_seq_id: str
    refget_sha512t24u: str
    refget_length: int
    # Populated during Phase B
    sequences_identical: bool | None = None
    seqrepo_seq_length: int | None = None
    refget_seq_length: int | None = None
    has_stop_codon: bool | None = None
    stop_codon_positions: list[int] = field(default_factory=list)
    # Digest computed from seqrepo's fetched bytes (no normalization)
    seqrepo_raw_digest: str | None = None
    # Digest computed from refgetstore's fetched bytes (no normalization)
    refget_raw_digest: str | None = None
    # Whether seqrepo's stored seq_id matches seq_seqhash (normalized, * stripped)
    seqrepo_digest_matches_normalized: bool | None = None
    # Whether refgetstore's stored digest matches sha512t24u of its own bytes
    refget_digest_self_consistent: bool | None = None
    # Classification
    cause: str = ""  # "star_normalization" or "sequence_divergence"
    # Error during fetch
    fetch_error: str | None = None


@dataclass
class DiagnosisReport:
    """Overall report of the investigation."""

    total_ensembl_aliases_checked: int = 0
    aliases_both_resolve: int = 0
    aliases_digests_agree: int = 0
    aliases_digests_mismatch: int = 0
    aliases_only_in_refget: int = 0
    aliases_only_in_seqrepo: int = 0
    mismatches: list[MismatchRecord] = field(default_factory=list)
    # Phase C classification
    all_mismatches_are_star_normalization: bool | None = None
    root_cause: str = ""
    recommendation: str = ""


def get_seqrepo_seq_ids(sr: SeqRepo) -> dict[str, str]:
    """Map ensembl aliases to seqrepo seq_ids.

    Returns {alias: seq_id} for all current Ensembl aliases.
    """
    conn = sr.aliases._db  # noqa: SLF001
    cursor = conn.execute(
        "SELECT alias, seq_id FROM seqalias "
        "WHERE namespace = 'Ensembl' AND is_current = 1"
    )
    result: dict[str, str] = {}
    for alias, seq_id in cursor:
        result[alias] = seq_id
    return result


def get_refget_ensembl_digests(store_path: Path) -> dict[str, str]:
    """Map ensembl aliases to refgetstore sha512t24u digests.

    Reads the ensembl.tsv alias sidecar directly for bulk access.
    Format is ``alias\\tdigest`` per line.

    Returns {alias: sha512t24u} for all ensembl namespace aliases.
    """
    tsv_path = store_path / "aliases" / "sequences" / "ensembl.tsv"
    result: dict[str, str] = {}
    with tsv_path.open("r") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) == 2:
                result[parts[0]] = parts[1]
    return result


def phase_a(
    sr: SeqRepo, store: RefgetStore, store_path: Path, report: DiagnosisReport
) -> list[MismatchRecord]:
    """Phase A: Enumerate all mismatches between seqrepo and refgetstore."""
    print("\n=== Phase A: Enumerating ensembl alias mismatches ===")

    sr_aliases = get_seqrepo_seq_ids(sr)
    rg_aliases = get_refget_ensembl_digests(store_path)

    sr_keys = set(sr_aliases)
    rg_keys = set(rg_aliases)
    common = sr_keys & rg_keys

    report.total_ensembl_aliases_checked = len(sr_keys | rg_keys)
    report.aliases_only_in_seqrepo = len(sr_keys - rg_keys)
    report.aliases_only_in_refget = len(rg_keys - sr_keys)
    report.aliases_both_resolve = len(common)

    mismatches: list[MismatchRecord] = []
    agree_count = 0
    for alias in sorted(common):
        sr_id = sr_aliases[alias]
        rg_id = rg_aliases[alias]
        if sr_id == rg_id:
            agree_count += 1
        else:
            # Resolve length from refgetstore metadata
            rec = store.get_sequence_by_alias("ensembl", alias)
            length = rec.metadata.length if rec else -1
            mismatches.append(
                MismatchRecord(
                    accession=alias,
                    seqrepo_seq_id=sr_id,
                    refget_sha512t24u=rg_id,
                    refget_length=length,
                )
            )

    report.aliases_digests_agree = agree_count
    report.aliases_digests_mismatch = len(mismatches)
    report.mismatches = mismatches

    print(f"  Total aliases checked:   {report.total_ensembl_aliases_checked}")
    print(f"  Both resolve, agree:     {agree_count}")
    print(f"  Both resolve, MISMATCH:  {len(mismatches)}")
    print(f"  Only in seqrepo:         {report.aliases_only_in_seqrepo}")
    print(f"  Only in refgetstore:     {report.aliases_only_in_refget}")

    # Breakdown by prefix
    enst_mm = sum(1 for m in mismatches if m.accession.startswith("ENST"))
    ensp_mm = sum(1 for m in mismatches if m.accession.startswith("ENSP"))
    print(f"  ENST mismatches:         {enst_mm}")
    print(f"  ENSP mismatches:         {ensp_mm}")

    return mismatches


def phase_b(
    sr: SeqRepo, store: RefgetStore, mismatches: list[MismatchRecord]
) -> None:
    """Phase B: Fetch full sequences from both stores and compare byte-by-byte."""
    print(f"\n=== Phase B: Byte-level comparison of {len(mismatches)} mismatches ===")

    for i, mm in enumerate(mismatches):
        # Fetch from seqrepo
        try:
            sr_seq: str = sr.fetch_uri(f"Ensembl:{mm.accession}")
        except Exception as exc:  # noqa: BLE001
            mm.fetch_error = f"seqrepo fetch failed: {exc}"
            print(f"  [{i+1}/{len(mismatches)}] {mm.accession}: {mm.fetch_error}")
            continue

        # Fetch from refgetstore
        try:
            store.load_sequence(mm.refget_sha512t24u)
            rg_seq: str = store.get_substring(
                mm.refget_sha512t24u, 0, mm.refget_length
            )
        except Exception as exc:  # noqa: BLE001
            mm.fetch_error = f"refgetstore fetch failed: {exc}"
            print(f"  [{i+1}/{len(mismatches)}] {mm.accession}: {mm.fetch_error}")
            continue

        mm.seqrepo_seq_length = len(sr_seq)
        mm.refget_seq_length = len(rg_seq)
        mm.sequences_identical = sr_seq == rg_seq

        # Check for stop codon *
        mm.has_stop_codon = "*" in sr_seq or "*" in rg_seq
        mm.stop_codon_positions = [j for j, c in enumerate(sr_seq) if c == "*"]

        # Compute sha512t24u from each side's raw bytes (no normalization)
        mm.seqrepo_raw_digest = sha512t24u(sr_seq.encode("ascii"))
        mm.refget_raw_digest = sha512t24u(rg_seq.encode("ascii"))

        # Check if seqrepo's stored seq_id matches seq_seqhash (which strips *)
        normalized_digest: str = seq_seqhash(sr_seq, normalize=True)
        mm.seqrepo_digest_matches_normalized = mm.seqrepo_seq_id == normalized_digest

        # Check if refgetstore's stored digest is self-consistent with its own bytes
        mm.refget_digest_self_consistent = (
            mm.refget_sha512t24u == mm.refget_raw_digest
        )

        # Classify the cause
        if mm.sequences_identical and mm.has_stop_codon:
            mm.cause = "star_normalization"
        elif not mm.sequences_identical:
            mm.cause = "sequence_divergence"
        else:
            mm.cause = "unknown"

        if (i + 1) % 20 == 0 or (i + 1) == len(mismatches):
            print(f"  Processed {i+1}/{len(mismatches)}")


def phase_c(report: DiagnosisReport) -> None:
    """Phase C: Classify root cause across all mismatches."""
    print("\n=== Phase C: Root cause classification ===")

    mismatches = report.mismatches
    if not mismatches:
        print("  No mismatches to classify.")
        return

    # Count by cause
    by_cause: dict[str, list[MismatchRecord]] = {}
    for m in mismatches:
        by_cause.setdefault(m.cause, []).append(m)
    n_fetch_errors = sum(1 for m in mismatches if m.fetch_error is not None)

    total = len(mismatches)
    n_star = len(by_cause.get("star_normalization", []))
    n_diverge = len(by_cause.get("sequence_divergence", []))
    n_unknown = len(by_cause.get("unknown", []))
    n_sr_self = sum(1 for m in mismatches if m.seqrepo_digest_matches_normalized)
    n_rg_self = sum(1 for m in mismatches if m.refget_digest_self_consistent)

    print(f"  Total mismatches:                       {total}")
    print(f"  Cause: * normalization (same bytes):    {n_star}")
    print(f"  Cause: sequence divergence (diff bytes):{n_diverge}")
    print(f"  Cause: unknown:                         {n_unknown}")
    print(f"  Fetch errors:                           {n_fetch_errors}")
    print(f"  seqrepo id matches normalized digest:   {n_sr_self}/{total}")
    print(f"  refget id self-consistent:              {n_rg_self}/{total}")

    # Show divergence details
    if n_diverge > 0:
        print(f"\n  Sequence divergence details:")
        for m in by_cause["sequence_divergence"]:
            print(
                f"    {m.accession}: "
                f"seqrepo={m.seqrepo_seq_length}bp refget={m.refget_seq_length}bp "
                f"(different content across Ensembl releases)"
            )

    report.all_mismatches_are_star_normalization = n_star == total

    report.root_cause = (
        f"{n_star}/{total} mismatches are caused by "
        "bioutils.sequences.normalize_sequence() stripping '*' (stop codon) "
        "characters before computing the sha512t24u digest. SeqRepo uses this "
        "normalized digest as its seq_id (and thus its ga4gh:SQ.* alias), while "
        "gtars RefgetStore computes sha512t24u directly from the raw FASTA bytes "
        "including '*'. The underlying sequence bytes are identical in both stores "
        "for these cases."
    )
    if n_diverge > 0:
        diverged = [m.accession for m in by_cause["sequence_divergence"]]
        report.root_cause += (
            f" The remaining {n_diverge} mismatch(es) ({', '.join(diverged)}) "
            "are genuine sequence content differences between seqrepo's Ensembl "
            "load (from an older release) and refgetstore's Ensembl r113."
        )

    report.recommendation = (
        "The refgetstore digest is correct per the GA4GH VRS specification, "
        "which defines the sequence digest as sha512t24u over the raw sequence "
        "bytes without normalization. SeqRepo's digest for *-containing proteins "
        "is technically non-conformant. Document the divergence and list the "
        "affected accessions. For the sequence-divergent transcript(s), the "
        "refgetstore carries the r113 version which is authoritative for that release."
    )

    print(f"\n  Root cause: * normalization ({n_star}) + sequence divergence ({n_diverge})")


def phase_d(report: DiagnosisReport, report_path: Path) -> None:
    """Phase D: Write summary and JSON report."""
    print(f"\n=== Phase D: Summary ===")
    print(f"  Root cause: {report.root_cause[:120]}...")
    print(f"  Recommendation: {report.recommendation[:120]}...")

    # Show a few examples
    print("\n  Sample mismatches:")
    for mm in report.mismatches[:8]:
        stars = len(mm.stop_codon_positions)
        print(
            f"    {mm.accession:30s}  "
            f"len={mm.refget_length:>5d}  "
            f"stars={stars:>2d}  "
            f"cause={mm.cause}"
        )

    # Write JSON report
    output = asdict(report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n  Report written to {report_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refget-store",
        type=Path,
        default=DEFAULT_STORE_PATH,
        help="Path to refgetstore directory (default: %(default)s)",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=HERE / "diagnosis_report.json",
        help="Path for JSON report output (default: %(default)s)",
    )
    args = parser.parse_args()

    seqrepo_root = os.environ.get("SEQREPO_ROOT_DIR")
    if not seqrepo_root:
        print("ERROR: SEQREPO_ROOT_DIR not set", file=sys.stderr)
        return 1
    if not args.refget_store.exists():
        print(
            f"ERROR: refget store not found at {args.refget_store}",
            file=sys.stderr,
        )
        return 1

    print(f"seqrepo:      {seqrepo_root}")
    print(f"refgetstore:  {args.refget_store}")

    sr = SeqRepo(seqrepo_root)
    store = RefgetStore.open_local(str(args.refget_store))
    store.pull_aliases()

    report = DiagnosisReport()

    t0 = time.monotonic()
    mismatches = phase_a(sr, store, args.refget_store, report)
    phase_b(sr, store, mismatches)
    phase_c(report)
    phase_d(report, args.report)
    elapsed = time.monotonic() - t0

    print(f"\nCompleted in {elapsed:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
