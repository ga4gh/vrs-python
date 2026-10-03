#!/usr/bin/env python
"""Ensembl-specific translator equivalence: seqrepo vs refgetstore.

Verifies that VRS Allele IDs agree between backends for SPDI expressions
using Ensembl transcript accessions. This complements the bulk_equivalence.py
test (which uses genomic VCF inputs) by exercising the Ensembl alias
resolution code paths in the data proxy.

SPDI format is used instead of HGVS to avoid the UTA database dependency.
The SPDI path exercises: ``derive_refget_accession`` → ``get_sequence`` →
``normalize`` — the full alias-to-allele pipeline.

Usage:
    SEQREPO_ROOT_DIR=/path/to/seqrepo/2024-12-20 \\
        uv run python misc/refgetstore/ensembl_translator_equivalence.py \\
            [--refget-store misc/refgetstore/store]

Exits 0 if all Allele IDs match, 1 if any mismatch or error.
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path

from biocommons.seqrepo import SeqRepo

from ga4gh.vrs.dataproxy import CompositeDataProxy, DataProxy
from ga4gh.vrs.dataproxy_refget import (
    RefgetAliasDataProxy,
    RefgetSequenceDataProxy,
)
from ga4gh.vrs.dataproxy_seqrepo import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator

HERE = Path(__file__).resolve().parent
DEFAULT_STORE_PATH = HERE / "store"


@dataclass
class SPDICase:
    """A single SPDI test case with Ensembl accession."""

    spdi: str
    description: str


# Well-known, stable Ensembl transcripts present in both seqrepo 2024-12-20
# and refgetstore (r113). Versions were selected to exist in both stores.
# None of these contain '*' in their sequence (DNA transcripts), so the
# digest normalization difference does not apply.
SPDI_CASES: list[SPDICase] = [
    # BRCA1 — well-known cancer gene, both stores have v9
    SPDICase("ENST00000357654.9:100:1:T", "BRCA1 SNV at 100"),
    SPDICase("ENST00000357654.9:5000:1:C", "BRCA1 SNV near 3' end"),
    SPDICase("ENST00000357654.9:200:0:A", "BRCA1 insertion at 200"),
    SPDICase("ENST00000357654.9:300:3:G", "BRCA1 delins at 300"),
    # EGFR — both stores have v9
    SPDICase("ENST00000261584.9:500:1:A", "EGFR SNV at 500"),
    SPDICase("ENST00000261584.9:2000:2:CC", "EGFR MNV at 2000"),
    # KRAS — both stores have v11
    SPDICase("ENST00000288602.11:100:1:T", "KRAS SNV at 100"),
    # BRAF — both stores have v7
    SPDICase("ENST00000275493.7:1000:1:A", "BRAF SNV at 1000"),
    # RB1 — both stores have v8
    SPDICase("ENST00000378567.8:500:1:G", "RB1 SNV at 500"),
    # Test deletion (empty insertion)
    SPDICase("ENST00000357654.9:400:1:", "BRCA1 1bp deletion at 400"),
]


def make_seqrepo_proxy(seqrepo_root: Path) -> SeqRepoDataProxy:
    sr = SeqRepo(root_dir=str(seqrepo_root), fd_cache_size=None)
    return SeqRepoDataProxy(sr)


def make_refget_proxy(store_path: Path) -> CompositeDataProxy:
    from gtars.refget import RefgetStore

    store = RefgetStore.open_local(str(store_path))
    return CompositeDataProxy(
        RefgetSequenceDataProxy(store),
        RefgetAliasDataProxy(store),
    )


def translate_spdi(
    dp: DataProxy, spdi_expr: str, *, do_normalize: bool = True
) -> tuple[str | None, str | None]:
    """Translate an SPDI expression and return (allele_id, error)."""
    try:
        tlr = AlleleTranslator(dp, identify=True)
        allele = tlr.translate_from(spdi_expr, fmt="spdi")
        if allele is None:
            return None, "translate_from returned None"
        return allele.id, None
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"


class Report:
    def __init__(self) -> None:
        self.passed = 0
        self.failed = 0

    def check(self, name: str, ok: bool, details: str = "") -> bool:
        mark = "PASS" if ok else "FAIL"
        line = f"  [{mark}] {name}"
        if details:
            line += f" — {details}"
        print(line)
        if ok:
            self.passed += 1
        else:
            self.failed += 1
        return ok


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refget-store",
        type=Path,
        default=DEFAULT_STORE_PATH,
        help="Path to refgetstore directory (default: %(default)s)",
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

    sr_dp = make_seqrepo_proxy(Path(seqrepo_root))
    rg_dp = make_refget_proxy(args.refget_store)

    r = Report()

    # Phase 1: Verify both backends can resolve the Ensembl accessions
    print("\n=== Phase 1: Accession resolution ===")
    accessions_seen: set[str] = set()
    for case in SPDI_CASES:
        ac = case.spdi.split(":")[0]
        if ac in accessions_seen:
            continue
        accessions_seen.add(ac)

        sr_acc = sr_dp.derive_refget_accession(ac)
        rg_acc = rg_dp.derive_refget_accession(ac)
        r.check(
            f"derive_refget_accession({ac})",
            sr_acc is not None and rg_acc is not None and sr_acc == rg_acc,
            f"seqrepo={sr_acc} refget={rg_acc}",
        )

    # Phase 2: Translate SPDI expressions and compare Allele IDs
    print("\n=== Phase 2: SPDI → VRS Allele equivalence ===")
    for case in SPDI_CASES:
        sr_id, sr_err = translate_spdi(sr_dp, case.spdi)
        rg_id, rg_err = translate_spdi(rg_dp, case.spdi)

        if sr_err:
            r.check(f"{case.spdi} ({case.description})", False, f"seqrepo error: {sr_err}")
            continue
        if rg_err:
            r.check(f"{case.spdi} ({case.description})", False, f"refget error: {rg_err}")
            continue

        r.check(
            f"{case.spdi} ({case.description})",
            sr_id == rg_id,
            f"seqrepo={sr_id} refget={rg_id}" if sr_id != rg_id else f"id={sr_id}",
        )

    # Phase 3: Verify sequence content agreement for test accessions
    print("\n=== Phase 3: Sequence content spot-checks ===")
    for ac in sorted(accessions_seen):
        try:
            sr_len = sr_dp.get_sequence_length(ac)
            rg_len = rg_dp.get_sequence_length(ac)
        except KeyError:
            r.check(f"get_sequence_length({ac})", False, "KeyError on one or both backends")
            continue
        r.check(
            f"get_sequence_length({ac})",
            sr_len == rg_len,
            f"seqrepo={sr_len} refget={rg_len}",
        )
        # Compare first 100 bytes
        sr_head = sr_dp.get_sequence(ac, start=0, end=min(100, sr_len))
        rg_head = rg_dp.get_sequence(ac, start=0, end=min(100, rg_len))
        r.check(
            f"get_sequence({ac})[0:100]",
            sr_head == rg_head,
            f"{'match' if sr_head == rg_head else 'MISMATCH'}",
        )

    print()
    print(f"{r.passed} passed, {r.failed} failed")
    print("OVERALL:", "EQUIVALENT" if r.failed == 0 else "DIVERGENT")
    return 0 if r.failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
