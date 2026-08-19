#!/usr/bin/env python
"""Verify the custom refgetstore built by build_store.py.

Checks that every sequence referenced by every assembly report was loaded,
that the alias graph contains only the namespaces vrs-python needs, and
that a handful of digests match the seqrepo 2024-12-20 ground truth.

Exits 0 on full success, 1 on any failure. Each check prints PASS/FAIL
on its own line; final line is OVERALL PASS or OVERALL FAIL.
"""

from __future__ import annotations

import sys
from pathlib import Path

from gtars.refget import RefgetStore

HERE = Path(__file__).resolve().parent
STORE_PATH = HERE / "store"
DOWNLOAD_DIR = HERE / "downloads"

EXPECTED_SEQ_NAMESPACES = {
    "GRCh38",
    "GRCh38.p14",
    "GRCh37",
    "GRCh37.p13",
    "refseq",
    "insdc",
    "ensembl",
}
EXPECTED_COLL_NAMESPACES = {"refseq", "insdc"}
FORBIDDEN_SEQ_NAMESPACES = {
    "sha512t24u",
    "ga4gh",
    "SEGUID",
    "SHA1",
    "VMC",
    "NCBI",
    "MD5",
}

# (namespace, alias) -> expected sha512t24u, captured from seqrepo 2024-12-20
# and the refgenie jungle store during earlier analysis in this session.
GROUND_TRUTH: dict[tuple[str, str], str] = {
    ("refseq", "NC_000001.11"): "Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
    ("refseq", "NC_000013.11"): "_0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
    ("GRCh38", "chr1"): "Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
    ("GRCh38", "chrM"): "k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct",
    ("GRCh38", "chrX"): "w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
    ("GRCh38", "chrY"): "8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
    ("GRCh38", "chr22"): "7B7SHsmchAR0dFcDCuSFjJAo7tX87krQ",
    ("GRCh38.p14", "chr1"): "Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
    ("insdc", "CM000663.2"): "Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
}

EXPECTED_COLL_REFSEQ_ALIASES = {
    "GCF_000001405.26",
    "GCF_000001405.40",
    "GCF_000001405.13",
    "GCF_000001405.25",
}

# (namespace, alias) -> (sha512t24u, length). Captured from the build_store.py
# run that ingested Ensembl release 113 FASTAs. These are the canonical
# latest versions in release 113 for a few well-known genes across all three
# biotypes we load (cdna, ncrna, pep). Bumping the Ensembl release will
# likely invalidate some of these — update after the rebuild.
ENSEMBL_GROUND_TRUTH: dict[tuple[str, str], tuple[str, int]] = {
    # cdna (protein-coding transcript)
    ("ensembl", "ENST00000256474.3"): ("xBKOKptLLDr-k4hTyCetvARn16pDS_rW", 4414),  # VHL
    ("ensembl", "ENST00000357654.9"): ("oR9jzdHf6J23TozeyuvTXtwlm6PpUHHl", 7088),  # BRCA1
    # ncrna (non-coding transcript)
    ("ensembl", "ENST00000429829.6"): ("FHnP3_NSzX0WZLjoHIPFjv89fSl18Xwk", 19245),  # XIST lncRNA
    # pep (protein)
    ("ensembl", "ENSP00000256474.3"): ("z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu", 213),   # VHL
    ("ensembl", "ENSP00000350283.3"): ("nUzIPnHMyQV52hzgBbKl5vlbSwx8M8_Y", 1863),  # BRCA1
    ("ensembl", "ENSP00000269305.4"): ("KAxM06sYzBF6zFftFaYq9E_18wsnn7al", 393),   # TP53
}


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


def parse_assembly_report(
    path: Path,
) -> tuple[str | None, str | None, list[dict[str, str]]]:
    """Mini copy of build_store.parse_assembly_report kept inline so this
    script is self-contained.
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
    return refseq_accn, genbank_accn, rows


def section_inventory(store: RefgetStore, r: Report) -> None:
    print("\n[A] Store inventory")
    stats = store.stats()
    print(f"    stats: {stats}")
    # 4 assembly collections + 16 RNA shards + 16 protein shards = 36. Lower
    # bound of 4 keeps the assembly-only configuration passing too.
    r.check(
        "n_collections >= 4",
        int(stats["n_collections"]) >= 4,
        f"got {stats['n_collections']}",
    )
    n_seq = int(stats["n_sequences"])
    # 941 observed on the first full build (Phase 1). Heavy dedup across
    # GRCh38/p14 and GRCh37/p13 keeps this much lower than the naive sum.
    r.check(
        "n_sequences >= 900",
        n_seq >= 900,
        f"got {n_seq}",
    )
    seq_ns = set(store.list_sequence_alias_namespaces())
    r.check(
        "sequence namespaces == expected set",
        seq_ns == EXPECTED_SEQ_NAMESPACES,
        f"got {sorted(seq_ns)}",
    )
    coll_ns = set(store.list_collection_alias_namespaces())
    r.check(
        "collection namespaces == expected set",
        coll_ns == EXPECTED_COLL_NAMESPACES,
        f"got {sorted(coll_ns)}",
    )


def section_forbidden(store: RefgetStore, r: Report) -> None:
    print("\n[B] Forbidden namespace check")
    seq_ns = set(store.list_sequence_alias_namespaces())
    intersection = FORBIDDEN_SEQ_NAMESPACES & seq_ns
    r.check(
        "no forbidden sequence alias namespaces",
        len(intersection) == 0,
        f"intersection: {sorted(intersection)}",
    )


def section_ground_truth(store: RefgetStore, r: Report) -> None:
    print("\n[C] Ground-truth digest cross-check (vs seqrepo 2024-12-20)")
    for (ns, alias), expected in GROUND_TRUTH.items():
        rec = store.get_sequence_by_alias(ns, alias)
        actual = rec.metadata.sha512t24u if rec is not None else None
        r.check(
            f"{ns}:{alias}",
            actual == expected,
            f"expected {expected}, got {actual}",
        )
    # metadata fidelity
    rec = store.get_sequence_by_alias("refseq", "NC_000001.11")
    if rec is None:
        r.check("NC_000001.11 metadata length", False, "record not found")
        r.check("NC_000001.11 metadata md5", False, "record not found")
    else:
        r.check(
            "refseq:NC_000001.11 length == 248956422",
            rec.metadata.length == 248956422,
            f"got {rec.metadata.length}",
        )
        r.check(
            "refseq:NC_000001.11 md5 matches",
            rec.metadata.md5 == "6aef897c3d6ff0c78aff06ac189178dd",
            f"got {rec.metadata.md5}",
        )


def section_coverage(store: RefgetStore, r: Report) -> None:
    print("\n[D] Assembly report coverage")
    for namespace in ("GRCh38", "GRCh38.p14", "GRCh37", "GRCh37.p13"):
        report_path = DOWNLOAD_DIR / f"{namespace}.assembly_report.txt"
        if not report_path.exists():
            r.check(f"{namespace} report present", False, f"{report_path} missing")
            continue
        _, _, rows = parse_assembly_report(report_path)
        miss_refseq: list[str] = []
        miss_ucsc: list[str] = []
        miss_insdc_ns: list[str] = []
        miss_insdc_global: list[str] = []
        total = 0
        for row in rows:
            refseq_ac = row.get("RefSeq-Accn", "").strip()
            ucsc = row.get("UCSC-style-name", "").strip()
            genbank = row.get("GenBank-Accn", "").strip()
            # Scope mirrors build_store.add_aliases_for_row: we only alias
            # sequences that are actually in the RefSeq dataset, so skip
            # rows where RefSeq-Accn is "na" (GenBank-only contigs that
            # don't appear in the RefSeq _genomic.fna.gz).
            if not refseq_ac or refseq_ac.lower() == "na":
                continue
            total += 1
            if store.get_sequence_by_alias("refseq", refseq_ac) is None:
                miss_refseq.append(refseq_ac)
            if ucsc and ucsc.lower() != "na":
                if store.get_sequence_by_alias(namespace, ucsc) is None:
                    miss_ucsc.append(ucsc)
            if genbank and genbank.lower() != "na":
                if store.get_sequence_by_alias(namespace, genbank) is None:
                    miss_insdc_ns.append(genbank)
                if store.get_sequence_by_alias("insdc", genbank) is None:
                    miss_insdc_global.append(genbank)
        r.check(
            f"{namespace}: refseq coverage ({total} rows)",
            not miss_refseq,
            f"{len(miss_refseq)} missing, first: {miss_refseq[:3]}",
        )
        r.check(
            f"{namespace}: UCSC coverage",
            not miss_ucsc,
            f"{len(miss_ucsc)} missing, first: {miss_ucsc[:3]}",
        )
        r.check(
            f"{namespace}: GenBank coverage ({namespace}:)",
            not miss_insdc_ns,
            f"{len(miss_insdc_ns)} missing, first: {miss_insdc_ns[:3]}",
        )
        r.check(
            f"{namespace}: GenBank coverage (insdc:)",
            not miss_insdc_global,
            f"{len(miss_insdc_global)} missing, first: {miss_insdc_global[:3]}",
        )


def section_invariants(store: RefgetStore, r: Report) -> None:
    print("\n[E] Cross-patch / cross-build invariants")

    def digest_for(ns: str, alias: str) -> str | None:
        rec = store.get_sequence_by_alias(ns, alias)
        return rec.metadata.sha512t24u if rec is not None else None

    grch38_chr1 = digest_for("GRCh38", "chr1")
    p14_chr1 = digest_for("GRCh38.p14", "chr1")
    grch37_chr1 = digest_for("GRCh37", "chr1")
    p13_chr1 = digest_for("GRCh37.p13", "chr1")
    grch38_chrm = digest_for("GRCh38", "chrM")
    # GRCh37 base has no mitochondrion contig; p13 is the first to include it
    # as MT/chrM/NC_012920.1. Compare against p13 (and refseq:NC_012920.1).
    grch37p13_chrm = digest_for("GRCh37.p13", "chrM")
    refseq_mt = digest_for("refseq", "NC_012920.1")

    r.check(
        "GRCh38:chr1 == GRCh38.p14:chr1",
        grch38_chr1 is not None and grch38_chr1 == p14_chr1,
        f"{grch38_chr1} vs {p14_chr1}",
    )
    r.check(
        "GRCh37:chr1 == GRCh37.p13:chr1",
        grch37_chr1 is not None and grch37_chr1 == p13_chr1,
        f"{grch37_chr1} vs {p13_chr1}",
    )
    r.check(
        "GRCh38:chr1 != GRCh37:chr1",
        grch38_chr1 is not None
        and grch37_chr1 is not None
        and grch38_chr1 != grch37_chr1,
        f"{grch38_chr1} vs {grch37_chr1}",
    )
    r.check(
        "GRCh38:chrM == GRCh37.p13:chrM (mito stable across builds)",
        grch38_chrm is not None and grch38_chrm == grch37p13_chrm,
        f"{grch38_chrm} vs {grch37p13_chrm}",
    )
    r.check(
        "GRCh38:chrM == refseq:NC_012920.1 (mito accession round-trip)",
        grch38_chrm is not None and grch38_chrm == refseq_mt,
        f"{grch38_chrm} vs {refseq_mt}",
    )


def section_sequence_bytes(store: RefgetStore, r: Report) -> None:
    print("\n[F] Sequence data smoke test")
    rec = store.get_sequence_by_alias("refseq", "NC_000001.11")
    if rec is None:
        r.check("retrieve chr1 slice", False, "NC_000001.11 not found")
        return
    digest = rec.metadata.sha512t24u
    try:
        store.load_sequence(digest)  # type: ignore[attr-defined]
        data = store.get_substring(digest, 1_000_000, 1_000_050)
    except Exception as exc:  # noqa: BLE001
        r.check("retrieve chr1 slice", False, f"load/get_substring raised: {exc}")
        return
    r.check(
        "chr1 slice length == 50",
        data is not None and len(data) == 50,
        f"got len={len(data) if data else None}",
    )
    r.check(
        "chr1 slice is ACGTN-only",
        data is not None and set(data.upper()) <= set("ACGTN"),
        f"unexpected chars: {sorted(set(data.upper()) - set('ACGTN')) if data else None}",
    )


def section_collection_aliases(store: RefgetStore, r: Report) -> None:
    print("\n[G] Collection aliases")
    refseq = set(store.list_collection_aliases("refseq") or [])
    r.check(
        "refseq collection aliases >= expected set",
        EXPECTED_COLL_REFSEQ_ALIASES <= refseq,
        f"got {sorted(refseq)}",
    )
    insdc = store.list_collection_aliases("insdc") or []
    r.check(
        "insdc collection aliases has 4 entries",
        len(insdc) >= 4,
        f"got {sorted(insdc)}",
    )


def section_ensembl(store: RefgetStore, r: Report) -> None:
    """Ensembl coverage checks.

    Asserts that each canonical ENST/ENSP identifier in
    ``ENSEMBL_GROUND_TRUTH`` resolves, its metadata matches the captured
    digest and length, and that at least one sequence's bytes can be
    read via ``get_substring`` (smoke-tests the full read path through
    ``RefgetSequenceDataProxy._ensure_sequence_loaded``).
    """
    print("\n[H] Ensembl coverage smoke test")
    for (ns, alias), (expected_digest, expected_length) in ENSEMBL_GROUND_TRUTH.items():
        rec = store.get_sequence_by_alias(ns, alias)
        if rec is None:
            r.check(f"{ns}:{alias}", False, "record not found")
            continue
        md = rec.metadata
        r.check(
            f"{ns}:{alias} digest",
            md.sha512t24u == expected_digest,
            f"expected {expected_digest}, got {md.sha512t24u}",
        )
        r.check(
            f"{ns}:{alias} length == {expected_length}",
            md.length == expected_length,
            f"got {md.length}",
        )

    # Sequence bytes smoke test on VHL cDNA — proves the FASTA payload
    # made it to disk for an Ensembl-loaded seqset, not just the alias row.
    rec = store.get_sequence_by_alias("ensembl", "ENST00000256474.3")
    if rec is None:
        r.check("VHL cDNA byte fetch", False, "record not found")
        return
    digest = rec.metadata.sha512t24u
    try:
        store.load_sequence(digest)  # type: ignore[attr-defined]
        data = store.get_substring(digest, 0, 50)
    except Exception as exc:  # noqa: BLE001
        r.check("VHL cDNA byte fetch", False, f"load/get_substring raised: {exc}")
        return
    r.check(
        "VHL cDNA slice length == 50",
        data is not None and len(data) == 50,
        f"got len={len(data) if data else None}",
    )
    r.check(
        "VHL cDNA slice is ACGTN-only",
        data is not None and set(data.upper()) <= set("ACGTN"),
        f"unexpected chars: {sorted(set(data.upper()) - set('ACGTN')) if data else None}",
    )


KNOWN_DIVERGENT_PATH = HERE / "ensembl_known_divergent.txt"


def section_known_divergent(store: RefgetStore, r: Report) -> None:
    """Regression test for Ensembl accessions with known seqrepo digest divergence.

    Loads the list of accessions from ``ensembl_known_divergent.txt`` and
    verifies that each one still resolves in the store with the expected
    sha512t24u digest. A change here means the store was rebuilt with
    different Ensembl data and the divergent list needs refreshing.
    """
    print("\n[I] Known-divergent Ensembl accessions")
    if not KNOWN_DIVERGENT_PATH.exists():
        r.check("known-divergent file exists", False, f"not found at {KNOWN_DIVERGENT_PATH}")
        return

    entries: list[tuple[str, str]] = []
    with KNOWN_DIVERGENT_PATH.open("r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                # accession, seqrepo_seq_id, refget_sha512t24u, [cause]
                entries.append((parts[0], parts[2]))

    r.check(
        f"known-divergent entries loaded ({len(entries)})",
        len(entries) >= 100,
        f"got {len(entries)}",
    )

    # Spot-check first 10 entries: verify each resolves with expected digest
    checked = 0
    for accession, expected_digest in entries[:10]:
        rec = store.get_sequence_by_alias("ensembl", accession)
        if rec is None:
            r.check(f"{accession} resolves", False, "not found in store")
            continue
        r.check(
            f"{accession} digest unchanged",
            rec.metadata.sha512t24u == expected_digest,
            f"expected {expected_digest}, got {rec.metadata.sha512t24u}",
        )
        checked += 1

    r.check(
        f"spot-checked {checked} known-divergent entries",
        checked >= 10,
        f"checked {checked}",
    )


def main() -> int:
    if not STORE_PATH.exists():
        print(f"ERROR: store not found at {STORE_PATH}", file=sys.stderr)
        return 2
    store = RefgetStore.open_local(str(STORE_PATH))
    store.pull_aliases()

    r = Report()
    section_inventory(store, r)
    section_forbidden(store, r)
    section_ground_truth(store, r)
    section_coverage(store, r)
    section_invariants(store, r)
    section_sequence_bytes(store, r)
    section_collection_aliases(store, r)
    section_ensembl(store, r)
    section_known_divergent(store, r)

    print()
    print(f"{r.passed} passed, {r.failed} failed")
    print("OVERALL:", "PASS" if r.failed == 0 else "FAIL")
    return 0 if r.failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
