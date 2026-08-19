# misc/refgetstore

Build a custom `gtars.refget.RefgetStore` archive augmented with the
per-sequence and per-collection aliases vrs-python needs to use refget as
its alias backend (replacing seqrepo for alias mapping). This covers:

- Assembly aliases: `GRCh38:chr1`, `GRCh38.p14:chrM`, etc.
- Cross-assembly genomic accessions: `refseq:NC_000001.11`, `insdc:CM000663.2`.
- Collection-level: `refseq:GCF_*`, `insdc:GCA_*`.
- **NCBI RefSeq transcripts + proteins**: `refseq:NM_000551.3`,
  `refseq:NP_000542.1`, etc., loaded from NCBI RefSeq mRNA/Prot shards.
- **Ensembl transcripts + proteins**: `ensembl:ENST00000256474.3`,
  `ensembl:ENSP00000256474.3`, etc., loaded from Ensembl release cdna /
  ncrna / pep FASTAs (release is pinned in `assemblies.toml`).

No dependency on `bioutils` or `biocommons.seqrepo` — the loader reads NCBI
assembly reports directly and derives transcript/protein aliases from FASTA
header names.

## Prerequisites

- `gtars.refget` (already a vrs-python extras dependency).
- Network access to `ftp.ncbi.nlm.nih.gov` and `ftp.ensembl.org` (or
  pre-populated cache dirs).
- Disk:
  - ~5 GB for one GRCh38 major + GRCh38.p14 assemblies (the `.fna.gz` files
    are ~900 MB each; the Encoded store is ~2–3 GB).
  - ~3-5 GB additional for the RefSeq mRNA+protein shards (64 files,
    ~50-100 MB each).
  - ~500 MB additional for the Ensembl cdna + ncrna + pep FASTAs.

## Layout

    misc/refgetstore/
      assemblies.toml     # declarative assembly list (edit to add patches)
      build_store.py      # loader entry point
      README.md           # this file
      downloads/          # cached FASTA + assembly_report.txt (gitignored)
      store/              # output RefgetStore (gitignored)

## Usage

Build everything in `assemblies.toml`:

    uv run python misc/refgetstore/build_store.py

Build only one assembly:

    uv run python misc/refgetstore/build_store.py --assembly GRCh38

Build only one seqset:

    uv run python misc/refgetstore/build_store.py --seqset refseq_human_rna

Skip categories entirely:

    uv run python misc/refgetstore/build_store.py --skip-seqsets
    uv run python misc/refgetstore/build_store.py --skip-assemblies

All CLI flags:

    --config PATH         assemblies.toml (default: misc/refgetstore/assemblies.toml)
    --store-dir PATH      output RefgetStore (default: misc/refgetstore/store)
    --download-dir PATH   download cache (default: misc/refgetstore/downloads)
    --assembly NAME       only process this assembly namespace
    --seqset NAME         only process this seqset name
    --skip-assemblies     skip all [[assembly]] entries
    --skip-seqsets        skip all [[seqset]] entries
    --force-download      re-fetch even if cached files exist

## Config reference

Each `[[assembly]]` block:

| field | required | description |
| --- | --- | --- |
| `namespace` | yes | Sequence-alias namespace to populate (e.g. `GRCh38`, `GRCh38.p14`). Becomes the namespace of `GRCh38:chr1`-style aliases. |
| `fasta_url` | yes | URL of NCBI `*_genomic.fna.gz`. gtars ingests `.gz` directly. |
| `report_url` | yes | URL of the corresponding `*_assembly_report.txt`. |
| `load_fasta` | no, default `true` | If false, aliases are added only for digests already present in the store (cheap patch fanout). |
| `fasta_path` | no | Local path overriding the downloaded FASTA (relative to repo root). Used for the already-downloaded GRCh38 major. |

Each `[[seqset]]` block (flat FASTA where the header name is the accession):

| field | required | description |
| --- | --- | --- |
| `name` | yes | Human id for logs + CLI selection. |
| `namespace` | yes | Alias namespace for every ingested sequence (e.g. `refseq`, `ensembl`). |
| `url_template` | yes | URL with `{shard}` placeholder, or plain URL if unsharded. |
| `shard_range` | no | `[min, max]` inclusive substituted into `{shard}`. Omit for a single-file seqset. |

The same schema handles both the NCBI RefSeq mRNA/Prot shards (sharded
with a `{shard}` placeholder) and the Ensembl cdna / ncrna / pep
releases (single-file, `shard_range` omitted).

## What gets written

Per assembly run, the loader adds:

- One sequence collection per `load_fasta=true` assembly, named by its
  collection sha512t24u digest.
- Per-sequence aliases under `<namespace>` (the assembly namespace):
  - UCSC-style-name (`chr1`, `chrM`, …)
  - GenBank-Accn (`CM000663.2`, …)
- Per-sequence aliases under the cross-assembly `refseq` namespace
  (`NC_000001.11` → digest). Emitted once per digest.
- Per-sequence aliases under the cross-assembly `insdc` namespace
  (`CM000663.2` → digest). Emitted once per digest.
- Per-collection aliases under `refseq` (`GCF_000001405.40` → collection
  digest) and `insdc` (`GCA_000001405.29` → collection digest), parsed from
  the assembly report header.

Per seqset run, the loader adds:

- One sequence collection per FASTA file ingested (each
  `human.<shard>.rna.fna.gz` and each Ensembl `*.fa.gz` becomes its own
  collection).
- Per-sequence aliases under `<namespace>` using each FASTA header's first
  token. Examples:
  - NCBI RefSeq shards: `refseq:NM_000551.3`, `refseq:NP_000542.1`,
    `refseq:XM_*`, `refseq:XP_*`.
  - Ensembl cdna: `ensembl:ENST00000256474.3` (protein-coding transcripts).
  - Ensembl ncrna: `ensembl:ENST00000429829.6` (non-coding — XIST, miRNAs,
    snoRNAs, etc.).
  - Ensembl pep: `ensembl:ENSP00000256474.3`.

Aliases for `sha512t24u:…` and `ga4gh:SQ.…` are intentionally NOT written —
they are the raw digest with a prefix and are synthesized at query time by
`RefgetAliasDataProxy`.

## Version-drift policy

The refgetstore pins **one version** of each source:

- **NCBI assemblies**: pinned by GCF accession (e.g. `GCF_000001405.40` for
  GRCh38.p14). Each patch is a superset of the prior, so only the latest
  patch of each major assembly is loaded.
- **NCBI RefSeq transcripts/proteins**: the FTP shards
  (`human.*.rna.fna.gz` / `human.*.protein.faa.gz`) publish only the
  current version of each accession. Older versions are not available.
- **Ensembl**: pinned to a single release (currently **r113**, set in
  `assemblies.toml`). Only the latest version of each ENST/ENSP in that
  release is included.

This differs from seqrepo's accumulative model, which retains every
version it has ever loaded across multiple historical imports. Callers that
depend on older transcript/protein versions will encounter `KeyError` on
accessions the refgetstore has not loaded.

### Known digest divergence on Ensembl proteins with `*`

For ~0.1% of Ensembl protein sequences (115 ENSPs containing embedded
stop codon `*` characters), the refgetstore and seqrepo produce different
`ga4gh:SQ.*` identifiers **despite storing identical sequence bytes**.
This is caused by `bioutils.sequences.normalize_sequence()` stripping `*`
before computing the digest in seqrepo, while gtars computes `sha512t24u`
directly from the raw FASTA bytes. The refgetstore digest is correct per
the GA4GH VRS specification. See `ensembl_known_divergent.txt` and
`diagnosis_report.json` for the full list and investigation details.

DNA sequences are not affected — the divergence only applies to protein
sequences containing `*`.

## Verification scripts

| Script | Purpose |
| --- | --- |
| `verify_store.py` | 67+ post-build checks: inventory, ground-truth digests, coverage, Ensembl, known-divergent regression |
| `compare_against_seqrepo.py` | 3-phase functional equivalence against seqrepo 2024-12-20 |
| `vcf_equivalence.py` | VCF annotator backend comparison + timing |
| `bulk_equivalence.py` | Translator backend comparison on 1M VCF records |
| `ensembl_translator_equivalence.py` | SPDI-based translator comparison using Ensembl transcript accessions |
| `diagnose_ensembl_mismatches.py` | Root-cause analysis of the 116 Ensembl digest mismatches |
