#!/usr/bin/env python
"""Compare the refget-backed DataProxy against the seqrepo-backed DataProxy.

Verifies that the RefgetAliasDataProxy (backed by the custom refgetstore
archive under misc/refgetstore/store) behaves equivalently to the existing
SeqRepoDataProxy for the operations vrs-python callers actually use.

Three phases:

    Phase 1 — Functional equivalence on identifiers both proxies support.
        For every canonical identifier, seqrepo and refget must agree on:
            derive_refget_accession(x)
            translate_sequence_identifier(x, namespace="ga4gh")
            translate_sequence_identifier(x, namespace="refseq")
            get_sequence_length(x)

    Phase 2 — Capability extension. Identifiers where seqrepo raises
        KeyError but refget resolves cleanly (new namespaces / forms we
        added that seqrepo never had). Assert refget returns the expected
        digest.

    Phase 3 — Minimal get_aliases containment. get_aliases() returns full
        alias graphs that differ by design between the two backends
        (short-name vs UCSC, insdc, patch versions, legacy namespaces).
        For each canonical identifier, we check only that the
        refseq:<accn> and ga4gh:SQ.<digest> entries BOTH proxies emit
        are present on both sides — the only alias shapes vrs-python
        actually consumes via translate_sequence_identifier.

Usage:

    SEQREPO_ROOT_DIR=/path/to/seqrepo/2024-12-20 \\
        uv run python misc/refgetstore/compare_against_seqrepo.py

Exits 0 on full success, 1 on any mismatch.
"""

from __future__ import annotations

import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path

from ga4gh.vrs.dataproxy import DataProxy, create_dataproxy

HERE = Path(__file__).resolve().parent
REFGET_STORE_PATH = HERE / "store"


@dataclass
class CanonicalInput:
    identifier: str
    expected_refget_accession: str
    expected_refseq: str
    expected_length: int


# Inputs that BOTH proxies should resolve identically. Expected values come
# from earlier ground-truth spot-checks against seqrepo 2024-12-20.
CANONICAL_INPUTS: list[CanonicalInput] = [
    CanonicalInput("refseq:NC_000001.11", "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO", "NC_000001.11", 248956422),
    CanonicalInput("NC_000001.11",        "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO", "NC_000001.11", 248956422),
    CanonicalInput("refseq:NC_000013.11", "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT", "NC_000013.11", 114364328),
    CanonicalInput("refseq:NC_000023.11", "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP", "NC_000023.11", 156040895),
    CanonicalInput("refseq:NC_012920.1",  "SQ.k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct", "NC_012920.1",  16569),
    CanonicalInput("GRCh38:chr1",         "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO", "NC_000001.11", 248956422),
    CanonicalInput("GRCh38:chrM",         "SQ.k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct", "NC_012920.1",  16569),
    CanonicalInput("GRCh38:chr22",        "SQ.7B7SHsmchAR0dFcDCuSFjJAo7tX87krQ", "NC_000022.11", 50818468),
    CanonicalInput("GRCh38:chrX",         "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP", "NC_000023.11", 156040895),
    CanonicalInput("GRCh37:chr1",         "SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU", "NC_000001.10", 249250621),
    CanonicalInput("ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                   "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO", "NC_000001.11", 248956422),
]


@dataclass
class RefgetOnlyInput:
    identifier: str
    expected_refget_accession: str
    expected_length: int


# Inputs where seqrepo raises KeyError today but refget should resolve
# cleanly because we deliberately populated the relevant namespace.
#
# Bare ``SQ.<digest>`` is deliberately NOT tested here. VRS models store
# refgetAccession as ``SQ.<digest>`` internally, but every vrs-python code
# path that reads a refgetAccession and passes it to the dataproxy
# prepends ``ga4gh:`` first (see normalize.py:123, translator.py:453,
# hgvs_tools.py:223). Seqrepo has always been broken on bare ``SQ.<d>``
# inputs (``coerce_namespace`` raises ``ValueError``) and no caller
# constructs that form, so handling it is out of scope.
REFGET_ONLY_INPUTS: list[RefgetOnlyInput] = [
    # seqrepo 2024-12-20 does not carry a GRCh38.p14 namespace
    RefgetOnlyInput("GRCh38.p14:chr1", "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO", 248956422),
    # GRCh37.p13 IS in seqrepo (check GRCh37.p13:chr1) — this one is a sanity
    # check that our store also covers it.
    RefgetOnlyInput("GRCh37.p13:chr1", "SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU", 249250621),
    # GenBank cross-assembly namespace — seqrepo has no `insdc` namespace
    RefgetOnlyInput("insdc:CM000663.2", "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO", 248956422),
    # sha512t24u: prefix form — seqrepo also stores this as a real alias
    # row, but running it through refget alone keeps the phase test
    # simple and still exercises the _extract_digest prefix branch.
    RefgetOnlyInput("sha512t24u:Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO", 248956422),
]


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


def try_call(fn, *args, **kwargs):
    try:
        return fn(*args, **kwargs), None
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"


def phase1_equivalence(
    seq_dp: DataProxy, ref_dp: DataProxy, r: Report
) -> None:
    print("\n=== Phase 1: functional equivalence on shared inputs ===")
    for case in CANONICAL_INPUTS:
        ident = case.identifier
        print(f"\n-- {ident!r}")

        # derive_refget_accession: both must return the known digest.
        s_acc, _ = try_call(seq_dp.derive_refget_accession, ident)
        r_acc, _ = try_call(ref_dp.derive_refget_accession, ident)
        r.check(
            "derive_refget_accession agrees with expected",
            s_acc == case.expected_refget_accession
            and r_acc == case.expected_refget_accession,
            f"expected={case.expected_refget_accession!r} seqrepo={s_acc!r} refget={r_acc!r}",
        )

        # translate ns=ga4gh
        s_ga, _ = try_call(seq_dp.translate_sequence_identifier, ident, "ga4gh")
        r_ga, _ = try_call(ref_dp.translate_sequence_identifier, ident, "ga4gh")
        expected_ga = {f"ga4gh:{case.expected_refget_accession}"}
        r.check(
            "translate ns=ga4gh agrees with expected",
            set(s_ga or []) == expected_ga and set(r_ga or []) == expected_ga,
            f"seqrepo={sorted(s_ga or [])} refget={sorted(r_ga or [])}",
        )

        # translate ns=refseq
        s_rs, _ = try_call(seq_dp.translate_sequence_identifier, ident, "refseq")
        r_rs, _ = try_call(ref_dp.translate_sequence_identifier, ident, "refseq")
        expected_rs = {f"refseq:{case.expected_refseq}"}
        r.check(
            "translate ns=refseq agrees with expected",
            set(s_rs or []) == expected_rs and set(r_rs or []) == expected_rs,
            f"seqrepo={sorted(s_rs or [])} refget={sorted(r_rs or [])}",
        )

        # get_sequence_length
        s_len, _ = try_call(seq_dp.get_sequence_length, ident)
        r_len, _ = try_call(ref_dp.get_sequence_length, ident)
        r.check(
            "get_sequence_length agrees with expected",
            s_len == case.expected_length and r_len == case.expected_length,
            f"expected={case.expected_length} seqrepo={s_len} refget={r_len}",
        )


def phase2_refget_extensions(ref_dp: DataProxy, r: Report) -> None:
    print("\n=== Phase 2: refget-only capability extensions ===")
    for case in REFGET_ONLY_INPUTS:
        ident = case.identifier
        print(f"\n-- {ident!r}")
        r_acc, r_err = try_call(ref_dp.derive_refget_accession, ident)
        r.check(
            "derive_refget_accession",
            r_acc == case.expected_refget_accession,
            f"expected={case.expected_refget_accession!r} got={r_acc!r} err={r_err}",
        )
        r_len, r_err = try_call(ref_dp.get_sequence_length, ident)
        r.check(
            "get_sequence_length",
            r_len == case.expected_length,
            f"expected={case.expected_length} got={r_len} err={r_err}",
        )


def phase3_alias_containment(
    seq_dp: DataProxy, ref_dp: DataProxy, r: Report
) -> None:
    print("\n=== Phase 3: get_aliases containment for refseq: and ga4gh: ===")
    for case in CANONICAL_INPUTS:
        ident = case.identifier
        print(f"\n-- {ident!r}")
        s_all, _ = try_call(seq_dp.get_aliases, ident)
        r_all, _ = try_call(ref_dp.get_aliases, ident)
        s_set = set(s_all or [])
        r_set = set(r_all or [])
        expected_ga = f"ga4gh:{case.expected_refget_accession}"
        expected_rs = f"refseq:{case.expected_refseq}"
        r.check(
            "ga4gh:<digest> present on both sides",
            expected_ga in s_set and expected_ga in r_set,
            f"seqrepo={expected_ga in s_set} refget={expected_ga in r_set}",
        )
        r.check(
            "refseq:<accn> present on both sides",
            expected_rs in s_set and expected_rs in r_set,
            f"seqrepo={expected_rs in s_set} refget={expected_rs in r_set}",
        )


def main() -> int:
    # Suppress noisy DataProxy stack traces from seqrepo KeyErrors — they're
    # expected for phase 2 inputs and the comparison script handles them.
    logging.getLogger("ga4gh.vrs.dataproxy").setLevel(logging.CRITICAL)

    seqrepo_root = os.environ.get("SEQREPO_ROOT_DIR")
    if not seqrepo_root:
        print("ERROR: SEQREPO_ROOT_DIR not set", file=sys.stderr)
        return 2
    if not REFGET_STORE_PATH.exists():
        print(f"ERROR: refget store not found at {REFGET_STORE_PATH}", file=sys.stderr)
        return 2

    print(f"seqrepo: {seqrepo_root}")
    print(f"refget:  {REFGET_STORE_PATH}")

    seq_dp = create_dataproxy(f"seqrepo+file://{seqrepo_root}")
    ref_dp = create_dataproxy(f"refget+file://{REFGET_STORE_PATH}")

    r = Report()
    phase1_equivalence(seq_dp, ref_dp, r)
    phase2_refget_extensions(ref_dp, r)
    phase3_alias_containment(seq_dp, ref_dp, r)

    print()
    print(f"{r.passed} passed, {r.failed} failed")
    print("OVERALL:", "EQUIVALENT" if r.failed == 0 else "DIVERGENT")
    return 0 if r.failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
