"""Bulk translator-path equivalence test: seqrepo vs refgetstore.

Drives ga4gh.vrs.extras.translator.AlleleTranslator end-to-end over every
record of a large VCF with each dataproxy backend, then compares the
resulting VRS Allele IDs for exact agreement.

This complements ``vcf_equivalence.py`` (which drives the VCF *annotator*
and compares VCF output). ``bulk_equivalence.py`` drives the *translator*
directly via gnomad-style expressions so the test exercises the same
data-proxy call sites as HGVS/SPDI import paths do.

Seqrepo is always run with ``fd_cache_size=None`` (unbounded file-descriptor
cache) so the comparison reflects steady-state performance rather than
per-record file-open cost. Export a larger ``SEQREPO_LRU_CACHE_MAXSIZE``
(default 1_000_000) if your workload does more distinct sequence fetches
than that.

Everything runs single-threaded (no pysam thread pool, no VRS parallelism).

Usage:
    SEQREPO_ROOT_DIR=/path/to/seqrepo/2024-12-20 \\
      uv run python misc/refgetstore/bulk_equivalence.py \\
        --input ~/dev/vrs-rust/gnomad_chr1.1M.vcf.bgz \\
        --refget-store misc/refgetstore/store

The refget store must already be built by build_store.py.

Exits non-zero if any Allele ID mismatch is found. Writes a JSON summary
to ``--report`` if provided.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from biocommons.seqrepo import SeqRepo

from ga4gh.vrs.dataproxy import CompositeDataProxy
from ga4gh.vrs.dataproxy_refget import (
    RefgetAliasDataProxy,
    RefgetSequenceDataProxy,
)
from ga4gh.vrs.dataproxy_seqrepo import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator

if TYPE_CHECKING:
    from gtars.refget import RefgetStore


# ---------------------------------------------------------------------------
# Dataproxy construction (kept in sync with vcf_equivalence.py)
# ---------------------------------------------------------------------------


def open_refget_store(store_path: Path) -> RefgetStore:
    from gtars.refget import RefgetStore

    if not (store_path / "rgstore.json").exists():
        raise SystemExit(
            f"refget store not found at {store_path}; "
            f"run `uv run python misc/refgetstore/build_store.py` first"
        )
    return RefgetStore.open_local(str(store_path))


def make_seqrepo_proxy(seqrepo_root: Path) -> SeqRepoDataProxy:
    """Construct a seqrepo-backed proxy with an unbounded fd cache.

    ``SEQREPO_FD_CACHE_MAXSIZE`` env var overrides this if set; otherwise
    ``fd_cache_size=None`` means `functools.lru_cache(maxsize=None)` — no
    eviction at any point, so the seqrepo baseline measures steady state.
    """
    sr = SeqRepo(root_dir=str(seqrepo_root), fd_cache_size=None)
    return SeqRepoDataProxy(sr)


def make_refget_native_proxy(store: RefgetStore) -> CompositeDataProxy:
    return CompositeDataProxy(
        RefgetSequenceDataProxy(store),
        RefgetAliasDataProxy(store),
    )


# ---------------------------------------------------------------------------
# Translation + comparison
# ---------------------------------------------------------------------------


@dataclass
class BackendResult:
    name: str
    allele_ids: list[str | None] = field(default_factory=list)
    errors: list[str | None] = field(default_factory=list)
    wall_seconds: float = 0.0

    @property
    def n_success(self) -> int:
        return sum(1 for e in self.errors if e is None)

    @property
    def n_errors(self) -> int:
        return len(self.errors) - self.n_success


def iter_gnomad_expressions(input_vcf: Path, limit: int | None):
    """Yield (line_no, gnomad_expr, coord_str) for each record/alt-allele.

    One input VCF record with N alt alleles produces N gnomad expressions —
    matching how annotators split multi-allelic sites.
    """
    import pysam

    vcf = pysam.VariantFile(str(input_vcf))
    i = 0
    for rec in vcf:
        if not rec.alts:
            continue
        for alt in rec.alts:
            i += 1
            expr = f"{rec.chrom}-{rec.pos}-{rec.ref}-{alt}"
            coord = expr
            yield i, expr, coord
            if limit is not None and i >= limit:
                vcf.close()
                return
    vcf.close()


def run_backend(
    name: str,
    tlr: AlleleTranslator,
    input_vcf: Path,
    limit: int | None,
    progress_every: int,
) -> BackendResult:
    result = BackendResult(name=name)
    t0 = time.perf_counter()
    for i, expr, _coord in iter_gnomad_expressions(input_vcf, limit):
        try:
            allele = tlr.translate_from(expr, fmt="gnomad")
            result.allele_ids.append(allele.id if allele is not None else None)
            result.errors.append(None if allele is not None else "returned None")
        except Exception as exc:
            result.allele_ids.append(None)
            result.errors.append(f"{type(exc).__name__}: {exc}")
        if progress_every and i % progress_every == 0:
            elapsed = time.perf_counter() - t0
            rate = i / elapsed if elapsed > 0 else 0.0
            print(
                f"  [{name}] {i} records in {elapsed:.1f}s ({rate:.0f}/s)",
                file=sys.stderr,
            )
    result.wall_seconds = time.perf_counter() - t0
    return result


@dataclass
class Mismatch:
    index: int
    expr: str
    baseline_id: str | None
    other_id: str | None
    baseline_error: str | None
    other_error: str | None


def compare_results(
    baseline: BackendResult,
    other: BackendResult,
    exprs: list[str],
    max_mismatches: int,
) -> tuple[int, int, list[Mismatch]]:
    total = min(len(baseline.allele_ids), len(other.allele_ids))
    id_matches = 0
    mismatches: list[Mismatch] = []
    for i in range(total):
        b_id = baseline.allele_ids[i]
        o_id = other.allele_ids[i]
        if b_id == o_id and baseline.errors[i] == other.errors[i]:
            id_matches += 1
        elif len(mismatches) < max_mismatches:
            mismatches.append(
                Mismatch(
                    index=i,
                    expr=exprs[i],
                    baseline_id=b_id,
                    other_id=o_id,
                    baseline_error=baseline.errors[i],
                    other_error=other.errors[i],
                )
            )
    return total, id_matches, mismatches


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    home = Path.home()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=home / "dev/vrs-rust/gnomad_chr1.1M.vcf.bgz",
        help="Input VCF (bgz/plain). Default: %(default)s",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Stop after N alt-allele records (default: process entire file).",
    )
    parser.add_argument(
        "--seqrepo-root",
        type=Path,
        default=Path(os.environ.get("SEQREPO_ROOT_DIR", "")),
    )
    parser.add_argument(
        "--refget-store",
        type=Path,
        default=Path("misc/refgetstore/store"),
    )
    parser.add_argument(
        "--assembly",
        default="GRCh38",
        help="Default assembly name for the translator.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help="Optional JSON report path.",
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=20,
        help="Max mismatches to collect for reporting.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=5000,
        help="Emit a progress line every N records (0 to disable).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.input.exists():
        print(f"Error: input not found: {args.input}", file=sys.stderr)
        return 2
    if not args.seqrepo_root or not args.seqrepo_root.exists():
        print(
            f"Error: seqrepo root not found: {args.seqrepo_root}",
            file=sys.stderr,
        )
        return 2

    print("=" * 70)
    print("Bulk equivalence: translator path")
    print("=" * 70)
    print(f"  Input:        {args.input}")
    print(f"  SeqRepo root: {args.seqrepo_root}")
    print(f"  Refget store: {args.refget_store}")
    print(f"  Assembly:     {args.assembly}")
    print(f"  Limit:        {args.limit if args.limit else 'all records'}")
    print(
        f"  SEQREPO_FD_CACHE_MAXSIZE={os.environ.get('SEQREPO_FD_CACHE_MAXSIZE', '(unset, constructor=None=unbounded)')}"
    )
    print(
        f"  SEQREPO_LRU_CACHE_MAXSIZE={os.environ.get('SEQREPO_LRU_CACHE_MAXSIZE', '(unset, default=1_000_000)')}"
    )
    print()

    print("[setup] Opening refget store...")
    store = open_refget_store(args.refget_store)

    print("[setup] Building dataproxies...")
    seqrepo_proxy = make_seqrepo_proxy(args.seqrepo_root)
    refget_native_proxy = make_refget_native_proxy(store)

    seqrepo_tlr = AlleleTranslator(seqrepo_proxy, default_assembly_name=args.assembly)
    refget_native_tlr = AlleleTranslator(
        refget_native_proxy, default_assembly_name=args.assembly
    )

    # Collect expressions once for mismatch reporting
    print("[setup] Indexing expressions from input VCF...")
    exprs = [expr for _, expr, _ in iter_gnomad_expressions(args.input, args.limit)]
    print(f"         {len(exprs)} expressions to translate")
    print()

    backends: dict[str, BackendResult] = {}

    print("[run] seqrepo...")
    backends["seqrepo"] = run_backend(
        "seqrepo", seqrepo_tlr, args.input, args.limit, args.progress_every
    )
    print(f"      {backends['seqrepo'].wall_seconds:.2f}s")

    print("[run] refget_native...")
    backends["refget_native"] = run_backend(
        "refget_native",
        refget_native_tlr,
        args.input,
        args.limit,
        args.progress_every,
    )
    print(f"      {backends['refget_native'].wall_seconds:.2f}s")

    # --- Compare against seqrepo baseline ---
    baseline = backends["seqrepo"]
    print("\n" + "=" * 70)
    print("Results")
    print("=" * 70)
    print(
        f"  baseline (seqrepo):     successes={baseline.n_success} "
        f"errors={baseline.n_errors} wall={baseline.wall_seconds:.2f}s"
    )

    report_payload: dict = {
        "input": str(args.input),
        "limit": args.limit,
        "total_expressions": len(exprs),
        "backends": {},
    }

    all_match = True
    for name, result in backends.items():
        report_payload["backends"][name] = {
            "wall_seconds": result.wall_seconds,
            "n_success": result.n_success,
            "n_errors": result.n_errors,
        }
        if name == "seqrepo":
            continue
        total, id_matches, mismatches = compare_results(
            baseline, result, exprs, args.max_mismatches
        )
        print(
            f"  {name:28s}: id_matches={id_matches}/{total} "
            f"wall={result.wall_seconds:.2f}s "
            f"(ratio={baseline.wall_seconds / result.wall_seconds:.2f}x)"
        )
        report_payload["backends"][name]["id_matches"] = id_matches
        report_payload["backends"][name]["total_compared"] = total
        report_payload["backends"][name]["mismatches"] = [asdict(m) for m in mismatches]
        if id_matches != total:
            all_match = False
            n_show = min(5, len(mismatches))
            print(f"    First {n_show} mismatches:")
            for m in mismatches[:n_show]:
                print(
                    f"      [{m.index}] {m.expr}: "
                    f"{m.baseline_id} (err={m.baseline_error}) vs "
                    f"{m.other_id} (err={m.other_error})"
                )

    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report_payload, indent=2))
        print(f"\n  report written to {args.report}")

    if all_match:
        print("\nPASS: all backends agree on Allele IDs")
        return 0
    print("\nFAIL: at least one backend diverged from seqrepo baseline")
    return 1


if __name__ == "__main__":
    sys.exit(main())
