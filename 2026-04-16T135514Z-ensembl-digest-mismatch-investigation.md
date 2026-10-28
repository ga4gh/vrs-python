# Investigation: 116 Ensembl Digest Mismatches Between SeqRepo and RefgetStore

**Date**: 2026-04-16
**Branch**: `kf/dataproxy-backend-refactor`
**Investigator**: Kyle Ferriter + Claude

## 1. Background

The `kf/dataproxy-backend-refactor` branch introduces a refget-backed
data proxy as an alternative to seqrepo. A cross-check of the `ensembl`
namespace between seqrepo (2024-12-20) and the custom refgetstore
(Ensembl r113) found **116 aliases** where the stored digest identifiers
diverge despite both systems being loaded from Ensembl protein/transcript
FASTAs.

Totals from the alias-level cross-check:

```
total aliases checked:    600,587
both resolve, agree:      322,286
both resolve, MISMATCH:       116  (115 ENSP + 1 ENST)
only in seqrepo:           65,833
only in refgetstore:      212,352
```

This report documents the investigation, root cause, and impact.

## 2. Initial Hypotheses

Five hypotheses were ranked before investigation (from the
[comparison report](2026-04-15T211800Z-refgetstore-seqrepo-comparison.md)):

1. Ensembl silently re-published proteins under the same version number
2. **Seqrepo's loader stripped something (e.g. `*`) that refgetstore kept**
3. Cross-release contamination in seqrepo
4. Isoform equivalence class ambiguity
5. Bug in the ingest path

## 3. Investigation Method

### 3.1 Direct byte comparison

The first diagnostic step was to fetch the sequence bytes from both stores
for a sample mismatch (ENSP00000340011.8, KIR2DS4) and compare them:

```python
from biocommons.seqrepo import SeqRepo
from gtars.refget import RefgetStore
from ga4gh.core.digests import sha512t24u

sr = SeqRepo('/path/to/seqrepo/2024-12-20')
sr_seq = sr.fetch_uri('Ensembl:ENSP00000340011.8')

store = RefgetStore.open_local('misc/refgetstore/store')
rec = store.get_sequence_by_alias('ensembl', 'ENSP00000340011.8')
store.load_sequence(rec.metadata.sha512t24u)
rg_seq = store.get_substring(rec.metadata.sha512t24u, 0, rec.metadata.length)
```

**Result**: The sequences are **byte-identical** (297 aa, both uppercase,
both containing `*` at positions 239, 269, 280). Yet the stored
identifiers differ:

```
seqrepo seq_id:          fv4uQW4AfrwGj3dx9D4xfpccU7Qnqtia
refgetstore sha512t24u:  8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0
```

Computing `sha512t24u` from the raw bytes confirms only the refgetstore
digest:

```python
>>> sha512t24u(sr_seq.encode('ascii'))
'8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0'  # matches refgetstore, NOT seqrepo
```

### 3.2 Tracing the digest divergence

Two different functions compute what appear to be the same hash but
produce different results:

**GA4GH VRS digest** (`src/ga4gh/core/digests.py:12-34`):

```python
def sha512t24u(blob: bytes) -> str:
    digest_size = 24
    digest = hashlib.sha512(blob).digest()
    tdigest_b64us = base64.urlsafe_b64encode(digest[:digest_size])
    return tdigest_b64us.decode("ascii")
```

This takes raw bytes and hashes them directly. No normalization.

**bioutils digest** (`bioutils.digests.seq_seqhash`):

```python
def seq_seqhash(seq, normalize=True):
    seq = normalize_sequence(seq) if normalize else seq
    return str(vmc_digest(seq, digest_size=24))
```

This calls `normalize_sequence()` before hashing.

**The normalization function** (`bioutils.sequences.normalize_sequence`):

```python
def normalize_sequence(seq):
    """Removes whitespace and asterisks, and uppercases the string."""
    nseq = re.sub(r"[\s\*]", "", seq).upper()
    m = re.search("[^A-Z]", nseq)
    if m:
        raise RuntimeError("Normalized sequence contains non-alphabetic characters")
    return nseq
```

The key line: `re.sub(r"[\s\*]", "", seq)` — **strips all `*` characters**
before hashing.

### 3.3 How SeqRepo uses the normalized digest

SeqRepo computes its internal `seq_id` using `seq_seqhash()` (with
normalization enabled by default). The `seq_id` becomes:
- The primary key for sequence storage
- The value behind the `ga4gh:SQ.*` alias namespace

When `SeqRepoDataProxy.get_aliases()` (`dataproxy_seqrepo.py:51-59`)
returns aliases for a sequence, the `ga4gh:SQ.{seq_id}` entry uses this
normalized digest. Any caller that calls `derive_refget_accession()` gets
back `SQ.{normalized_digest}`.

### 3.4 How RefgetStore uses the raw digest

gtars `RefgetStore` computes `sha512t24u` directly from the FASTA bytes.
When `RefgetAliasDataProxy.get_aliases()` (`dataproxy_refget.py:216-251`)
synthesizes the `ga4gh:SQ.*` alias, it uses this raw digest:

```python
# dataproxy_refget.py:249-250
aliases.append(f"sha512t24u:{digest}")
aliases.append(f"ga4gh:SQ.{digest}")
```

### 3.5 Proof: DNA sequences are unaffected

For a DNA sequence like NC_000001.11 (chr1), `normalize_sequence()` has
no effect — DNA bases have no `*` to strip, and uppercasing is a no-op
for already-uppercase FASTA. Both algorithms produce identical digests:

```
NC_000001.11 seqrepo seq_id:  Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO
NC_000001.11 ga4gh sha512t24u: Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO  ← identical
```

This is why all 464,830 refseq namespace aliases and all 322,286
non-`*`-containing Ensembl aliases agree perfectly.

## 4. Automated Verification

A diagnostic script (`misc/refgetstore/diagnose_ensembl_mismatches.py`)
was written to verify the hypothesis across all 116 mismatches. It
operates in 4 phases:

1. **Phase A**: Enumerate all `ensembl` namespace aliases from both stores
   (reading seqrepo's SQLite directly and refgetstore's
   `aliases/sequences/ensembl.tsv`), identify where digests differ.

2. **Phase B**: For each mismatch, fetch full sequences from both stores,
   compare byte-by-byte, compute `sha512t24u` from each side's raw bytes,
   and check whether seqrepo's stored `seq_id` matches
   `seq_seqhash(seq, normalize=True)`.

3. **Phase C**: Classify each mismatch as `star_normalization` (identical
   bytes + `*` present) or `sequence_divergence` (different bytes).

4. **Phase D**: Output JSON report + human-readable summary.

### Results

```
Total mismatches:                        116
Cause: * normalization (same bytes):     115
Cause: sequence divergence (diff bytes):   1
Fetch errors:                              0
seqrepo id matches normalized digest:    116/116
refget id self-consistent:               116/116
```

**Every** seqrepo `seq_id` matches the result of `seq_seqhash(bytes)`.
**Every** refgetstore digest matches `sha512t24u(bytes)`. The two hash
functions simply disagree on what to hash (raw bytes vs `*`-stripped bytes).

### Stop codon distribution

The 115 `star_normalization` mismatches contain between 1 and 23 internal
stop codons:

| `*` count | Proteins |
|-----------|----------|
| 1         | 54       |
| 2         | 19       |
| 3         | 17       |
| 4         | 3        |
| 6-23      | 22       |

### Representative examples

```
ENSP00000375609.1  len= 225  stars= 1  pos=[222]
  seqrepo (normalized): 0wFdmvBugnimf4tOFjCVdckpw3P6GmOx
  refget  (raw):        Qg7FS7qyO-pfHp7HYd4wvr26DEx1oov4
  sha512t24u(bytes):    Qg7FS7qyO-pfHp7HYd4wvr26DEx1oov4  ← matches refget

ENSP00000340011.8  len= 297  stars= 3  pos=[239, 269, 280]
  seqrepo (normalized): fv4uQW4AfrwGj3dx9D4xfpccU7Qnqtia
  refget  (raw):        8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0
  sha512t24u(bytes):    8Co4fICavUFj9qH_GGXd5cbWbzh8oSZ0  ← matches refget

ENSP00000434417.1  len=2459  stars=17  pos=[2130, 2153, ..., 2418]
  seqrepo (normalized): zKGzcqTwIi1cpsh_25jTR-g7cycsIEpf
  refget  (raw):        Ih9FQHBcIpxLvdvSmQwf1B_uW07yCbjE
  sha512t24u(bytes):    Ih9FQHBcIpxLvdvSmQwf1B_uW07yCbjE  ← matches refget
```

In every case, `sha512t24u(actual_bytes)` matches the refgetstore digest,
not the seqrepo `seq_id`.

## 5. The Sequence-Divergent Outlier

ENST00000668831.1 is a genuine data difference, not a hashing issue:

```
seqrepo:      1,214 bp  (loaded 2020-04-12, from an older Ensembl release)
refgetstore:  1,120 bp  (Ensembl r113)
First base:   seqrepo='T' vs refget='A' — completely different content
```

Both stores are self-consistent (their stored digests match
`sha512t24u` of their own bytes). Ensembl re-annotated this transcript
between releases under the same version number — a rare but documented
occurrence for complex genomic loci.

## 6. Which Side Is Correct?

The GA4GH VRS specification defines the sequence identifier as
`sha512t24u` computed over the raw sequence bytes. The specification
does not mandate stripping `*` or any other normalization before
hashing. Therefore:

- **RefgetStore is correct** per GA4GH VRS: it computes `sha512t24u`
  directly from the stored bytes.
- **SeqRepo is non-conformant** for protein sequences containing `*`:
  its `seq_id` (and thus its `ga4gh:SQ.*` alias) is computed from
  `*`-stripped bytes, producing a different identifier.

However, this non-conformance is embedded in the `bioutils` library
(`bioutils.digests.seq_seqhash` → `bioutils.sequences.normalize_sequence`)
which is upstream of both seqrepo and vrs-python. Changing it would be a
breaking change affecting all existing seqrepo identifiers for
`*`-containing sequences.

## 7. Practical Impact

### 7.1 DNA workloads: zero impact

DNA sequences never contain `*`. All genomic operations (VCF annotation,
HGVS translation, allele normalization) are completely unaffected.
This was confirmed by:
- 1,000,000 VCF records: zero allele ID mismatches
  (`misc/refgetstore/bulk_equivalence.py`)
- 74/74 canonical identifier checks
  (`misc/refgetstore/compare_against_seqrepo.py`)

### 7.2 Ensembl transcript workloads: zero impact (for transcripts without `*`)

SPDI expressions using Ensembl transcript accessions (ENST*) produce
identical VRS Allele IDs on both backends. Confirmed by 25/25 checks
across BRCA1, EGFR, KRAS, BRAF, and RB1 transcripts
(`misc/refgetstore/ensembl_translator_equivalence.py`).

### 7.3 Protein sequences with `*`: identifier divergence

If a caller:
1. Resolves an ENSP accession containing `*` to a `ga4gh:SQ.*` digest
2. Uses that digest as a VRS identifier or passes it to another system

...they will get different identifiers depending on which backend was used.
This affects 115 out of 123,845 Ensembl protein entries (0.09%).

In practice, vrs-python's translator primarily works with genomic and
transcript-level variants, not protein sequences directly. The SPDI and
gnomAD translation paths use chromosome/transcript accessions, not ENSP.
HGVS protein expressions (`p.Val600Glu`) are resolved through UTA at
the transcript level, not via direct protein sequence lookup.

## 8. Artifacts Produced

| File | Purpose |
| --- | --- |
| `misc/refgetstore/diagnose_ensembl_mismatches.py` | 4-phase diagnostic script (Phase A-D) |
| `misc/refgetstore/diagnosis_report.json` | Full JSON report with all 116 mismatch records |
| `misc/refgetstore/ensembl_known_divergent.txt` | TSV of all 116 accessions with cause classification |
| `misc/refgetstore/ensembl_translator_equivalence.py` | SPDI-based Ensembl translator equivalence test |
| `misc/refgetstore/verify_store.py` section `[I]` | Regression test for known-divergent digests |
| `misc/refgetstore/README.md` | Version-drift policy and `*` divergence documentation |

## 9. Verification Summary

| Check | Result |
| --- | --- |
| `verify_store.py` (67 checks incl. new `[I]` section) | 67 pass, 0 fail |
| `compare_against_seqrepo.py` (74 canonical checks) | 74 pass, 0 fail |
| `test_dataproxy_refget.py` (29 unit tests) | 29 pass |
| `ensembl_translator_equivalence.py` (25 checks) | 25 pass, 0 fail |

## 10. Conclusion

The 116 "digest mismatches" break down as:

- **115/116 (99.1%)**: Not a data issue. Caused by
  `bioutils.sequences.normalize_sequence()` stripping `*` characters
  before `seq_seqhash()` computes the SHA-512 digest. The actual
  sequence bytes in both stores are identical. RefgetStore's digest is
  GA4GH VRS-conformant; seqrepo's is a legacy artifact of bioutils
  normalization. No fix needed — documented and regression-tested.

- **1/116 (0.9%)**: Genuine sequence content difference for
  ENST00000668831.1 across Ensembl releases (seqrepo loaded from ~2020,
  refgetstore from r113 2024). Both stores are self-consistent. The
  refgetstore version is authoritative for r113.

The refgetstore is a safe drop-in replacement for seqrepo across all
DNA/transcript workloads. The `*` normalization divergence is limited to
protein sequences and does not affect any known vrs-python code path.
