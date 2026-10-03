# Refgetstore vs seqrepo 2024-12-20 — coverage comparison snapshot

Snapshot date: 2026-04-15. Captured after the `ensembl_human_cdna`,
`ensembl_human_ncrna`, and `ensembl_human_pep` seqsets landed on top of
the existing 4 assemblies + 32 NCBI RefSeq mRNA/Prot shards.

This doc inventories what's in our custom refgetstore
(`misc/refgetstore/store/`), compares it alias-by-alias to the seqrepo
instance at `/Users/kferrite/dev/data/seqrepo/2024-12-20`, and flags
the one investigation thread that's still open (Ensembl ENSP digest
mismatches).

## 1. Refgetstore inventory

Raw stats (`RefgetStore.stats()`):

```
n_sequences:       783 625     unique sha512t24u digests across all sources
n_collections:      39          4 assemblies + 32 RefSeq shards + 3 Ensembl FASTAs
storage_mode:      Encoded     (on-disk, lazy sequence byte loading)
```

Sources loaded (from [assemblies.toml](misc/refgetstore/assemblies.toml)):

| Source                                                   | Collection count | File count | Fresh aliases written |
|----------------------------------------------------------|-----------------:|-----------:|----------------------:|
| NCBI GRCh38 major (`GCF_000001405.26`)                   | 1                | 1          | 1 820                 |
| NCBI GRCh38.p14 (`GCF_000001405.40`)                     | 1                | 1          | 1 918                 |
| NCBI GRCh37 major (`GCF_000001405.13`)                   | 1                | 1          |   344                 |
| NCBI GRCh37.p13 (`GCF_000001405.25`)                     | 1                | 1          |   694                 |
| NCBI RefSeq human mRNA (`human.1..16.rna.fna.gz`)        | 16               | 16         | 273 558               |
| NCBI RefSeq human protein (`human.1..16.protein.faa.gz`) | 16               | 16         | 197 939               |
| Ensembl r113 cdna (`Homo_sapiens.GRCh38.cdna.all.fa.gz`) | 1                | 1          | 207 133               |
| Ensembl r113 ncrna (`Homo_sapiens.GRCh38.ncrna.fa.gz`)   | 1                | 1          | 203 776               |
| Ensembl r113 pep (`Homo_sapiens.GRCh38.pep.all.fa.gz`)   | 1                | 1          | 123 845               |

Per-namespace alias counts (`list_sequence_alias_namespaces()`):

| Namespace       | Count   | Source                                        |
|-----------------|--------:|-----------------------------------------------|
| `GRCh38`        |     910 | GRCh38 major assembly report — chr*/CM*/etc.  |
| `GRCh38.p14`    |  1 410  | GRCh38.p14 assembly report                    |
| `GRCh37`        |     184 | GRCh37 major assembly report                  |
| `GRCh37.p13`    |     390 | GRCh37.p13 assembly report                    |
| `insdc`         |     941 | GenBank-Accn column across all assemblies     |
| `refseq`        | 472 438 | RefSeq-Accn across assemblies + NM/NR/XM/XR/NP/XP from mRNA/Prot shards |
| `ensembl`       | 534 754 | ENST/ENSP from Ensembl r113 cdna+ncrna+pep    |

Collection namespaces:

- `refseq:GCF_*` — 4 entries (the 4 NCBI assemblies)
- `insdc:GCA_*` — 4 entries (the 4 GenBank counterparts)

The seqset FASTAs (RefSeq shards and Ensembl files) each become their
own collection but do **not** get collection-level aliases — no
`refseq:GCF_*` equivalent exists for a transcript shard, and we don't
synthesize one.

## 2. Seqrepo inventory for comparison

Seqrepo namespaces at `/Users/kferrite/dev/data/seqrepo/2024-12-20`
(from `SELECT namespace, COUNT(*) FROM seqalias WHERE is_current=1 GROUP BY namespace`):

```
Ensembl                  388 235
NCBI                     885 113
GRCh37                       184
GRCh37.p2 ..  GRCh37.p13     256 .. 390
GRCh38                       910
GRCh38.p1 .. GRCh38.p12      926 .. 998
CHM1_1.1                     23
NCBI34                        3
NCBI35                       26
NCBI36                      149
JRGv1, JRGv2                 25  each
hs37-1kg                     84
hs37d5                       86
MD5, SEGUID, SHA1, VMC   1 144 093 each (digest-form aliases)
```

Things seqrepo has that our refgetstore does **not**:

- `GRCh38.p1` through `GRCh38.p13` — older patch levels (we pinned p14
  only since each patch is a superset of the prior).
- `GRCh37.p2`, `p5`, `p9`, `p10`, `p11`, `p12` — intermediate patches
  (we pinned p13 only).
- `CHM1_1.1`, `NCBI34`, `NCBI35`, `NCBI36`, `JRGv1`, `JRGv2`,
  `hs37-1kg`, `hs37d5` — legacy assemblies predating GRCh37/38, plus
  alt-contig sets.
- `MD5`, `SEGUID`, `SHA1`, `VMC` — digest-form aliases under non-SHA512t24u
  algorithms. These are derivable at query time, so we deliberately
  don't store them (see
  [build_store.py:add_aliases_for_row](misc/refgetstore/build_store.py)).

Things our refgetstore has that seqrepo does **not**:

- `GRCh38.p14` — our pinned patch level postdates what seqrepo has.
- `ensembl` **ncRNA biotype** (XIST, lincRNAs, miRNAs, snoRNAs,
  pseudogene transcripts, etc.) — seqrepo only loaded Ensembl cdna
  historically, so ~200 k ncRNA ENSTs exist only in our store.

## 3. Alias-level cross-check (refget vs seqrepo)

### 3.1 Assembly namespaces (genomic)

No systematic comparison was run; the canonical digests we care about
(chr1, chr22, chrX, chrY, chrM, NC_000001.11, CM000663.2) were already
verified by
[misc/refgetstore/verify_store.py](misc/refgetstore/verify_store.py)'s
sections [C] and [D]. 41/41 checks pass.

### 3.2 `refseq` / seqrepo `NCBI` namespace

Seqrepo calls this `NCBI`; our store calls it `refseq`. They contain
NM/NR/NP + NC/NG/NT/NW accessions.

| Bucket                                                 |    count | % of seqrepo NCBI |
|--------------------------------------------------------|---------:|------------------:|
| Both resolve, digests agree                            | 464 830  | 52.52 %           |
| Digest mismatch on the same versioned accession        |       0  |  0.00 %           |
| In seqrepo, not in refgetstore                         | 420 283  | 47.48 %           |
| In refgetstore, not in seqrepo                         |   7 608  |   n/a             |

**Zero digest mismatches**. The 420 k "seqrepo-only" rows are dominated
by historical versions that seqrepo has accumulated across loads but
that the current RefSeq FTP `mRNA_Prot/human.*.rna.fna.gz` /
`human.*.protein.faa.gz` shards no longer publish (NCBI ships "current"
only, the way Ensembl does). Our store matches whatever is in the
shards as of the download date.

Our `refseq` composition (by accession prefix):

```
XM_ : 127 028      # predicted mRNAs
XP_ : 127 028      # predicted proteins
NM_ :  70 898      # curated mRNAs
NP_ :  70 898      # curated proteins
XR_ :  51 132      # predicted ncRNAs
NR_ :  24 500      # curated ncRNAs
NW_ :     477      # genomic contigs from assemblies
NT_ :     415      # genomic contigs from assemblies
NC_ :      49      # chromosome/organelle entries from assemblies
YP_ :      13      # proteins from GenBank WGS, assembly-embedded
```

The `X*_` / `N*_` split matches the upstream RefSeq curation tiers
exactly: "model" (X*) vs "curated" (N*).

### 3.3 `ensembl` namespace

| Bucket                                                      |    count | % of seqrepo Ensembl |
|-------------------------------------------------------------|---------:|---------------------:|
| Both resolve, digests agree                                 | 322 286  | 83.01 %              |
| Digest mismatch on the same versioned accession (ENSP 115, ENST 1) |    116  |  0.03 %              |
| In seqrepo, not in refgetstore                              |  65 833  | 16.96 %              |
| In refgetstore, not in seqrepo                              | 212 352  |   n/a                |

83% exact overlap without any deliberate release alignment — that's
good confirmation that the seqset loader + our release-113 source is
sane. The 116 digest mismatches are discussed separately below.

Stem-level breakdown of the asymmetries (helpful for understanding
what's missing where):

**seqrepo-only (65 833 rows)** — versions or stems seqrepo has but
release 113 does not:

```
                              ENST     ENSP
same stem, diff version:     43 881    6 474     ← older transcript/protein
                                                   versions (BRCA1 .7 alongside .9)
stem not in refgetstore:      7 987    6 852     ← entire stems retired or never
                                                   in release 113
```

**refget-only (212 352 rows)** — versions or stems we have that
seqrepo never loaded:

```
                              ENST      ENSP
same stem, newer version:    22 862    4 732    ← normal bumps between seqrepo's
                                                  load date and release 113
stem not in seqrepo:        167 490   17 268    ← Ensembl ncRNA FASTA + new
                                                  genes/proteins added since 2024
```

The ~167 k "ENST stem not in seqrepo at all" rows are almost entirely
the `ncrna/` biotype, which seqrepo's historical loads omit — compare
with the 203 776 entries in `Homo_sapiens.GRCh38.ncrna.fa.gz`. About
30 k of those ncRNAs are also in seqrepo under some version (matching
the `same_stem` buckets); the rest (~167 k) are new coverage.

Breakdown of `ensembl` aliases by prefix:

| Prefix | refgetstore | seqrepo  |
|--------|------------:|---------:|
| `ENST` | 410 909     | 272 425  |
| `ENSP` | 123 845     | 115 171  |
| other  |           0 |     639  |

The ~639 "other" seqrepo Ensembl aliases are things like `CHR_*` and
`KI27*` — these are assembly contig aliases that somehow ended up
under the Ensembl namespace in an older seqrepo load. We don't
reproduce them.

## 4. Resolved — 116 Ensembl digest mismatches

**Status**: Investigated and root-caused (2026-04-16). See
[`misc/refgetstore/diagnose_ensembl_mismatches.py`](misc/refgetstore/diagnose_ensembl_mismatches.py)
and [`misc/refgetstore/diagnosis_report.json`](misc/refgetstore/diagnosis_report.json).

**Distribution**: 115 ENSP + 1 ENST.

### Root cause: two distinct mechanisms

#### 4a. `*` (stop codon) normalization — 115 ENSP mismatches

The underlying sequence bytes are **identical** in both stores.
The digest difference is caused by `bioutils.sequences.normalize_sequence()`
stripping `*` characters (and whitespace) before hashing:

- **SeqRepo** computes its `seq_id` via `bioutils.digests.seq_seqhash(seq)`,
  which calls `normalize_sequence(seq)` → strips `*` → uppercases →
  SHA-512 → truncate to 24 bytes → base64url. The resulting digest
  becomes the `ga4gh:SQ.{digest}` alias.
- **RefgetStore** (gtars) computes `sha512t24u` directly from the raw
  FASTA bytes, including `*`. This is the GA4GH VRS-conformant approach.

For protein sequences without `*`, the two algorithms produce the same
digest (as confirmed by the 322,286 agreeing aliases). Only the 115
ENSPs that contain embedded stop codons diverge.

Example (ENSP00000340011.8, KIR2DS4, 297 aa, 3 internal `*`):

```
seqrepo seq_id (normalized): fv4uQW4AfrwGj3dx9D4xfpccU7Qnqtia  ← sha512t24u("...TEQ TARILM...")
refget sha512t24u (raw):     8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0  ← sha512t24u("...TEQ*TARILM...")
                                                                    ^ this * is the difference
```

**Conclusion**: refgetstore's digest is correct per the GA4GH VRS
specification (sha512t24u of the raw bytes). SeqRepo's digest is a
legacy artifact of `bioutils.normalize_sequence()` stripping `*` before
hashing. This does NOT affect DNA sequences.

#### 4b. Sequence content divergence — 1 ENST mismatch

ENST00000668831.1 has genuinely different content: seqrepo loaded
1,214 bp (added 2020-04-12, from an older Ensembl release), while
refgetstore has 1,120 bp from release 113. The sequences share no
homology at position 0 — Ensembl re-annotated this transcript between
releases under the same version number.

Both stores are self-consistent (sha512t24u matches their own stored
bytes). The refgetstore carries the r113 version, which is
authoritative for that release.

### Impact on VRS interoperability

For vrs-python's primary use case (genomic DNA variants), this has
**zero impact** — DNA sequences never contain `*` and the digest
algorithms agree exactly. The divergence only affects callers who:

1. Look up protein sequences with embedded `*` (stop codons)
2. AND use the resulting `ga4gh:SQ.*` identifier for cross-system
   comparison between seqrepo and refgetstore

The 115 affected ENSPs are listed in
[`misc/refgetstore/ensembl_known_divergent.txt`](misc/refgetstore/ensembl_known_divergent.txt).

## 5. What this snapshot does NOT include

- No comparison of sequence bytes themselves — only digest equality
  was checked. If seqrepo and refgetstore both report the same
  sha512t24u, the bytes are identical by construction.
- No latency/throughput comparison beyond the earlier
  [bulk_equivalence.py](misc/refgetstore/bulk_equivalence.py) /
  [vcf_equivalence.py](misc/refgetstore/vcf_equivalence.py) runs,
  which showed refget ~1.15× faster than seqrepo on the translator
  path and ~1.14× on the VCF annotator path.
- No attempt to load GRCh37 Ensembl, Ensembl regulatory features
  (`ENSR*`), or Ensembl gene sequences (`ENSG*`). Those were
  explicitly out of scope for the release-113 seqset plan
  (`~/.claude/plans/rippling-splashing-cerf.md`).
- No check of seqrepo's non-current (`is_current=0`) rows. All counts
  use `WHERE is_current=1`.

## 6. Files of interest

- [misc/refgetstore/assemblies.toml](misc/refgetstore/assemblies.toml) —
  declarative source list, now with 4 assemblies + 2 sharded RefSeq
  seqsets + 3 Ensembl seqsets.
- [misc/refgetstore/build_store.py](misc/refgetstore/build_store.py) —
  loader; `process_seqset` uses the bulk `load_sequence_aliases` path
  that makes ~500 k new aliases finish in sub-second time.
- [misc/refgetstore/verify_store.py](misc/refgetstore/verify_store.py) —
  55 checks (was 41 before Ensembl) including a `[H]` Ensembl
  coverage section with VHL / BRCA1 / XIST / TP53 ground-truth.
- [misc/refgetstore/compare_against_seqrepo.py](misc/refgetstore/compare_against_seqrepo.py) —
  11 canonical-identifier checks (spot-check, not bulk); still 74/74
  passing after this round.
- `/Users/kferrite/dev/data/seqrepo/2024-12-20/aliases.sqlite3` — the
  seqrepo reference used for all cross-checks in this doc.

## 7. Suggested next steps (in rough priority order)

1. ~~**Investigate the 115 ENSP digest mismatches**~~ — **Done** (§4).
   Root cause: `bioutils.normalize_sequence()` strips `*` before
   hashing. RefgetStore digest is GA4GH-conformant; seqrepo's is not.
   No code fix needed — documented and accessions listed.
2. **Run a translator-path bulk_equivalence comparison that
   specifically stresses ENST/ENSP inputs** (e.g. feed
   `NM_*:c.1A>T`-style HGVS expressions that coerce into Ensembl via
   the translator's namespace fallback). The gnomAD chr1 VCF we
   currently use is purely genomic and doesn't exercise any of the
   Ensembl alias code paths.
3. **Decide whether to load GRCh37 Ensembl** (archive release 87).
   Probably still not needed, but worth revisiting if any caller
   reports KeyError on a GRCh37-only ENST.
4. **Document the version-drift reality in the README**: make clear
   that refgetstore ships "latest-per-release" only, not the
   cumulative version history that seqrepo accumulates. Callers that
   depend on older transcript versions need to load additional
   release snapshots.
5. **Consider a periodic rebuild script**. Pinning release 113 in
   `assemblies.toml` is good for reproducibility, but we should have
   a clear story for bumping it (update the URL, rerun, verify
   digests against ground-truth, update `ENSEMBL_GROUND_TRUTH` in
   `verify_store.py`).
