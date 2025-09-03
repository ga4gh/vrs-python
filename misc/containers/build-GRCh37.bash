#!/bin/bash

if [ -z $SEQREPO_ROOT_DIR ]; then
    echo "Must set SEQREPO_ROOT_DIR"
    exit 1
fi

# Load reference genome from pre-downloaded file
# File should already be present from Docker data layer
seqrepo -r $SEQREPO_ROOT_DIR init
seqrepo -r $SEQREPO_ROOT_DIR load -n NCBI GCF_000001405.13_GRCh37_genomic.fna.gz
seqrepo -r $SEQREPO_ROOT_DIR add-assembly-names