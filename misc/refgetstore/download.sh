#!/usr/bin/env bash
# Pre-fetch FASTAs and assembly reports into the loader's cache directory.
# Mirrors the filenames build_store.py expects so subsequent runs are cache hits.
# Supports resuming partial downloads via curl -C -.
set -euo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
DEST="${HERE}/downloads"
mkdir -p "${DEST}"

fetch() {
    local url="$1" out="$2"
    if [ -f "${out}" ]; then
        echo "[skip] ${out} (already present)"
        return
    fi
    echo "[get ] ${url}"
    curl -fL --retry 3 --retry-delay 5 -C - -o "${out}.part" "${url}"
    mv "${out}.part" "${out}"
}

base_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405"

# GRCh38 major (FASTA is already in the repo root via fasta_path override)
fetch "${base_url}/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_assembly_report.txt" \
      "${DEST}/GRCh38.assembly_report.txt"

# GRCh38.p14
fetch "${base_url}/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" \
      "${DEST}/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
fetch "${base_url}/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt" \
      "${DEST}/GRCh38.p14.assembly_report.txt"

# GRCh37 major
fetch "${base_url}/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.fna.gz" \
      "${DEST}/GCF_000001405.13_GRCh37_genomic.fna.gz"
fetch "${base_url}/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_assembly_report.txt" \
      "${DEST}/GRCh37.assembly_report.txt"

# GRCh37.p13
fetch "${base_url}/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz" \
      "${DEST}/GCF_000001405.25_GRCh37.p13_genomic.fna.gz"
fetch "${base_url}/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt" \
      "${DEST}/GRCh37.p13.assembly_report.txt"

echo
echo "Downloads ready in ${DEST}"
ls -lh "${DEST}"
