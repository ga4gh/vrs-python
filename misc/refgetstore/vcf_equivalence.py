"""Compare VRS annotation outputs from seqrepo vs refgetstore dataproxies.

Runs the vrs-python VCF annotator over the same input VCF with two
dataproxy configurations and reports match rates + end-to-end wall time:

  1. SeqRepoDataProxy (unbounded fd cache) — seqrepo baseline
  2. CompositeDataProxy(Refget, Refget-aliases) — fully refget-native

Seqrepo is always run with ``fd_cache_size=None`` (unbounded file-descriptor
cache) so the comparison reflects seqrepo's steady-state performance rather
than per-record open/close cost. ``SEQREPO_LRU_CACHE_MAXSIZE`` defaults to
1_000_000 in seqrepo; export a larger value if your input has more than
that many distinct sequence fetches so the sequence cache doesn't evict.

Everything runs single-threaded: pysam.VariantFile defaults to threads=1,
the translator path is Python-serial, and neither refget nor seqrepo uses
workers internally.

Usage:
    SEQREPO_ROOT_DIR=/path/to/seqrepo/2024-12-20 \\
      uv run python misc/refgetstore/vcf_equivalence.py \\
        --input-vcf ~/dev/vrs-rust/gnomad_chr1.1M.vcf.bgz \\
        --refget-store misc/refgetstore/store

The refget store must already be built by build_store.py. The script does
NOT ingest FASTA — it only opens existing stores.

Exits non-zero if any backend disagrees with the seqrepo baseline on
VRS_Allele_IDs.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from biocommons.seqrepo import SeqRepo

from ga4gh.vrs.dataproxy import CompositeDataProxy, DataProxy
from ga4gh.vrs.dataproxy_refget import (
    RefgetAliasDataProxy,
    RefgetSequenceDataProxy,
)
from ga4gh.vrs.dataproxy_seqrepo import SeqRepoDataProxy
from ga4gh.vrs.extras.annotator.vcf import VcfAnnotator

if TYPE_CHECKING:
    from gtars.refget import RefgetStore


# ---------------------------------------------------------------------------
# Dataproxy construction
# ---------------------------------------------------------------------------


def open_refget_store(store_path: Path) -> RefgetStore:
    """Open an existing on-disk refget store. Errors if not present."""
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
    """Composite: refget for both sequence bytes and alias lookups."""
    return CompositeDataProxy(
        RefgetSequenceDataProxy(store),
        RefgetAliasDataProxy(store),
    )


# ---------------------------------------------------------------------------
# Annotation + comparison
# ---------------------------------------------------------------------------


def annotate_vcf(
    dp: DataProxy,
    input_vcf: Path,
    output_vcf: Path,
    assembly: str,
) -> float:
    """Run VCF annotation and return wall-clock seconds."""
    annotator = VcfAnnotator(dp)
    t0 = time.perf_counter()
    annotator.annotate(
        input_vcf_path=input_vcf,
        output_vcf_path=output_vcf,
        vrs_attributes=True,
        assembly=assembly,
        compute_for_ref=True,
        require_validation=False,
    )
    return time.perf_counter() - t0


@dataclass
class CompareResult:
    total_records: int
    id_matches: int
    full_info_matches: int
    id_mismatches: list[str]


def compare_vcf_outputs(path_a: Path, path_b: Path) -> CompareResult:
    """Compare two annotated VCF files line-by-line on VRS INFO fields."""
    import pysam

    vcf_a = pysam.VariantFile(str(path_a))
    vcf_b = pysam.VariantFile(str(path_b))

    total = 0
    id_matches = 0
    id_mismatches: list[str] = []
    full_info_matches = 0

    vrs_fields = [
        "VRS_Allele_IDs",
        "VRS_Starts",
        "VRS_Ends",
        "VRS_States",
        "VRS_Lengths",
        "VRS_Error",
    ]

    for rec_a, rec_b in zip(vcf_a, vcf_b, strict=False):
        total += 1
        ids_a = rec_a.info.get("VRS_Allele_IDs", None)
        ids_b = rec_b.info.get("VRS_Allele_IDs", None)

        if ids_a == ids_b:
            id_matches += 1
        else:
            coord = (
                f"{rec_a.chrom}:{rec_a.pos}:{rec_a.ref}>{','.join(rec_a.alts or [])}"
            )
            id_mismatches.append(f"  {coord}: {ids_a} != {ids_b}")

        all_match = all(
            rec_a.info.get(f, None) == rec_b.info.get(f, None) for f in vrs_fields
        )
        if all_match:
            full_info_matches += 1

    vcf_a.close()
    vcf_b.close()

    return CompareResult(
        total_records=total,
        id_matches=id_matches,
        full_info_matches=full_info_matches,
        id_mismatches=id_mismatches,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    home = Path.home()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-vcf",
        type=Path,
        default=home / "dev/vrs-rust/gnomad_chr1.1M.vcf.bgz",
        help="Input VCF (bgz/plain). Default: %(default)s",
    )
    parser.add_argument(
        "--seqrepo-root",
        type=Path,
        default=Path(os.environ.get("SEQREPO_ROOT_DIR", "")),
        help="SeqRepo root dir (default: $SEQREPO_ROOT_DIR).",
    )
    parser.add_argument(
        "--refget-store",
        type=Path,
        default=Path("misc/refgetstore/store"),
        help="On-disk refget store built by build_store.py.",
    )
    parser.add_argument(
        "--assembly",
        default="GRCh38",
        help="Assembly passed to VcfAnnotator.annotate().",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(".pytest_tmp/vcf_equivalence"),
        help="Directory for annotated VCF outputs.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.input_vcf.exists():
        print(f"Error: input VCF not found: {args.input_vcf}", file=sys.stderr)
        return 2
    if not args.seqrepo_root or not args.seqrepo_root.exists():
        print(
            f"Error: seqrepo root not found: {args.seqrepo_root} "
            f"(set --seqrepo-root or $SEQREPO_ROOT_DIR)",
            file=sys.stderr,
        )
        return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_seqrepo = args.output_dir / "annotated_seqrepo.vcf"
    out_refget = args.output_dir / "annotated_refget.vcf"

    # --- Setup ---
    print("=" * 70)
    print("VCF Equivalence: seqrepo vs refgetstore")
    print("=" * 70)
    print(f"  Input VCF:    {args.input_vcf}")
    print(f"  SeqRepo root: {args.seqrepo_root}")
    print(f"  Refget store: {args.refget_store}")
    print(f"  Assembly:     {args.assembly}")
    print(f"  Output dir:   {args.output_dir}")
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

    timings: dict[str, float] = {}
    outputs: dict[str, Path] = {}

    # --- Annotations (single-threaded; pysam.VariantFile defaults to threads=1) ---

    print(f"\n[run] seqrepo (unbounded fd cache) -> {out_seqrepo}")
    t = annotate_vcf(seqrepo_proxy, args.input_vcf, out_seqrepo, args.assembly)
    timings["seqrepo"] = t
    outputs["seqrepo"] = out_seqrepo
    print(f"       done in {t:.2f}s")

    print(f"\n[run] refget native                -> {out_refget}")
    t = annotate_vcf(refget_native_proxy, args.input_vcf, out_refget, args.assembly)
    timings["refget_native"] = t
    outputs["refget_native"] = out_refget
    print(f"       done in {t:.2f}s")

    # --- Compare ---
    baseline_key = "seqrepo"
    baseline_path = outputs[baseline_key]

    print("\n" + "=" * 70)
    print("Results")
    print("=" * 70)

    all_match = True
    for key, path in outputs.items():
        if key == baseline_key:
            continue
        r = compare_vcf_outputs(baseline_path, path)
        total = r.total_records
        print(
            f"  {key:25s}  VRS_Allele_IDs={r.id_matches}/{total}  "
            f"full={r.full_info_matches}/{total}  "
            f"mismatches={len(r.id_mismatches)}"
        )
        if r.id_mismatches:
            all_match = False
            n_show = min(5, len(r.id_mismatches))
            for line in r.id_mismatches[:n_show]:
                print(line)
            if len(r.id_mismatches) > n_show:
                print(f"  ... and {len(r.id_mismatches) - n_show} more")

    print("\n  Timing:")
    base = timings[baseline_key]
    for key, t in timings.items():
        ratio = f"  ({base / t:.2f}x vs seqrepo)" if key != baseline_key else ""
        print(f"    {key:25s} {t:>8.2f}s{ratio}")

    if all_match:
        print("\nPASS: all backends agree on VRS_Allele_IDs")
        return 0
    print("\nFAIL: at least one backend diverged from seqrepo baseline")
    return 1


if __name__ == "__main__":
    sys.exit(main())
