#!/bin/bash
set -xeuo pipefail

if [ -z $SEQREPO_ROOT_DIR ]; then
    echo "Must set SEQREPO_ROOT_DIR"
    exit 1
fi

reference_url=https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.fna.gz
reference_fname=$(basename $reference_url)

download_reference() {
    curl -O $reference_url
    echo "$reference_url"
}

build_seqrepo() {
    # Load reference genome from pre-downloaded file
    # File should already be present from Docker data layer
    seqrepo -r $SEQREPO_ROOT_DIR init
    seqrepo -r $SEQREPO_ROOT_DIR load -n NCBI $reference_fname
    seqrepo -r $SEQREPO_ROOT_DIR add-assembly-names
}
