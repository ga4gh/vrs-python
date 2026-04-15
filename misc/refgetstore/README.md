# misc/refgetstore

Build a custom `gtars.refget.RefgetStore` archive augmented with the
per-sequence and per-collection aliases vrs-python needs to use refget as
its alias backend (replacing seqrepo for alias mapping). This covers:

- Assembly aliases: `GRCh38:chr1`, `GRCh38.p14:chrM`, etc.
- Cross-assembly genomic accessions: `refseq:NC_000001.11`, `insdc:CM000663.2`.
- Collection-level: `refseq:GCF_*`, `insdc:GCA_*`.
- **Transcripts + proteins**: `refseq:NM_000551.3`, `refseq:NP_000542.1`, etc.,
  loaded from NCBI RefSeq mRNA/Prot shards.

No dependency on `bioutils` or `biocommons.seqrepo` — the loader reads NCBI
assembly reports directly and derives transcript/protein aliases from FASTA
header names.

## Prerequisites

- `gtars.refget` (already a vrs-python extras dependency).
- Network access to `ftp.ncbi.nlm.nih.gov` (or pre-populated cache dirs).
- Disk:
  - ~5 GB for one GRCh38 major + GRCh38.p14 assemblies (the `.fna.gz` files
    are ~900 MB each; the Encoded store is ~2–3 GB).
  - ~3-5 GB additional for the RefSeq mRNA+protein shards (64 files,
    ~50-100 MB each).

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
| `namespace` | yes | Alias namespace for every ingested sequence (usually `refseq`). |
| `url_template` | yes | URL with `{shard}` placeholder, or plain URL if unsharded. |
| `shard_range` | no | `[min, max]` inclusive substituted into `{shard}`. Omit for a single-file seqset. |

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

- One sequence collection per shard FASTA (each `human.<shard>.rna.fna.gz`
  becomes its own collection).
- Per-sequence aliases under `<namespace>` using each FASTA header's first
  token (e.g. `refseq:NM_000551.3` → digest).

Aliases for `sha512t24u:…` and `ga4gh:SQ.…` are intentionally NOT written —
they are the raw digest with a prefix and are synthesized at query time by
`RefgetAliasDataProxy`.
