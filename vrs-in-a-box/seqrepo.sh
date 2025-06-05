#!/bin/bash
set -e  # Exit on any error


# Check if an argument is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fasta_file> <assembly_name>"
    exit 1
fi

FASTA_FILE=$1
ASSEMBLY_NAME=$2

export SEQREPO_ROOT=seqrepo-$ASSEMBLY_NAME

# Check if the file exists
if [ ! -f "$FASTA_FILE" ]; then
    echo "ERROR: File '$FASTA_FILE' not found!" 1>&2
    exit 1
fi

# Check if openrsync is in the output of /usr/bin/rsync --version
RSYNC_EXE=${RSYNC:-rsync}
if $RSYNC_EXE --version | grep -q "openrsync"; then
    if [ -x "/opt/homebrew/bin/rsync" ]; then
        RSYNC_OPTION="--rsync-exe /opt/homebrew/bin/rsync"
    else
        echo "ERROR: seqrepo requires 'rsync' (not 'openrsync') to run." 1>&2
        echo "       On MacOS rsync can be installed with Homebrew." 1>&2
        echo "       Set the RSYNC environment variable to specify the rsync executable to use." 1>&2
        exit 1
    fi
else
    RSYNC_OPTION="--rsync-exe $RSYNC_EXE"
fi

# Initialize seqrepo
echo "Initializing seqrepo in $SEQREPO_ROOT..."
seqrepo -r seqrepo --root-directory $SEQREPO_ROOT $RSYNC_OPTION init

# Load the provided FASTA file
echo "Loading FASTA file: $FASTA_FILE"
seqrepo -r seqrepo --root-directory $SEQREPO_ROOT $RSYNC_OPTION load "$FASTA_FILE" -n NCBI

# Add assembly names
echo "Adding assembly names..."
seqrepo -r seqrepo --root-directory $SEQREPO_ROOT $RSYNC_OPTION add-assembly-names

echo "Initial seqrepo build completed"

# vrs-annotate checks that GRCh38:1 is in the database as a sanity check
#   so for non-GRCh38, we have to add a dummy record so it will run
if [ "$ASSEMBLY_NAME" != "GRCh38" ]; then
sqlite3 $SEQREPO_ROOT/master/aliases.sqlite3 <<EOF
    INSERT INTO seqalias (seq_id, namespace, alias, added, is_current)
    VALUES ('Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO', 'GRCh38', '1', DATE('now'), 1);
EOF
sqlite3 $SEQREPO_ROOT/master/sequences/db.sqlite3 <<EOF
INSERT INTO seqinfo (seq_id, len, alpha, added, relpath)
VALUES ('Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO', 248956422, 'ACGMNRT', date('now'), '2025/0331/1411/1743430287.165611.fa.bgz');
EOF
echo "Spiked in GRCh38:1 sequence alias to make vrs-annotate run"
fi
# Create a tar.gz archive of the seqrepo directory
echo "Creating seqrepo archive..."
zip -r "${SEQREPO_ROOT}.zip" ${SEQREPO_ROOT}
rm -fR $SEQREPO_ROOT

echo "Archive created ${SEQREPO_ROOT}.zip"
