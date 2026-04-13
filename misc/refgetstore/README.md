# misc/refgetstore

Build a custom `gtars.refget.RefgetStore` archive augmented with the
per-sequence and per-collection aliases vrs-python needs to use refget as
its alias backend (replacing seqrepo for alias mapping). This covers the
shortlist: `GRCh38:chr1`, `refseq:NC_000001.11`, `insdc:CM000663.2`, plus the
collection-level `refseq:GCF_*` / `insdc:GCA_*` entries the jungle store
already publishes.

No dependency on `bioutils` or `biocommons.seqrepo` — the loader reads NCBI
assembly reports directly.

## Prerequisites

- `gtars.refget` (already a vrs-python extras dependency).
- Network access to `ftp.ncbi.nlm.nih.gov` (or pre-populated cache dirs).
- Disk: ~5 GB for one GRCh38 major + GRCh38.p14 (the .fna.gz files are ~900 MB
  each; the Encoded store is ~2–3 GB).

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

Build only one entry:

    uv run python misc/refgetstore/build_store.py --assembly GRCh38

All CLI flags:

    --config PATH        assemblies.toml (default: misc/refgetstore/assemblies.toml)
    --store-dir PATH     output RefgetStore (default: misc/refgetstore/store)
    --download-dir PATH  download cache (default: misc/refgetstore/downloads)
    --assembly NAME      only process this namespace from the config
    --force-download     re-fetch even if cached files exist

## Config reference

Each `[[assembly]]` block:

| field | required | description |
| --- | --- | --- |
| `namespace` | yes | Sequence-alias namespace to populate (e.g. `GRCh38`, `GRCh38.p14`). Becomes the namespace of `GRCh38:chr1`-style aliases. |
| `fasta_url` | yes | URL of NCBI `*_genomic.fna.gz`. gtars ingests `.gz` directly. |
| `report_url` | yes | URL of the corresponding `*_assembly_report.txt`. |
| `load_fasta` | no, default `true` | If false, aliases are added only for digests already present in the store (cheap patch fanout). |
| `fasta_path` | no | Local path overriding the downloaded FASTA (relative to repo root). Used for the already-downloaded GRCh38 major. |

## What gets written

Per run, the loader adds to the store:

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

Aliases for `sha512t24u:…` and `ga4gh:SQ.…` are intentionally NOT written —
they are the raw digest with a prefix and will be synthesized at query time
by the eventual `RefgetAliasDataProxy`.
