###
# podman build --arch linux/amd64,linux/arm64 --build-arg ASSEMBLY=GRCh38 -t docker.io/ga4gh/vrs-python:GRCh38 -f ./Dockerfile .
# podman build --arch linux/arm64 --build-arg ASSEMBLY=GRCh38 --target build -t docker.io/ga4gh/vrs-python:GRCh38-build -f ./Dockerfile .
###
# Data layer - downloads genomic reference files
FROM python:3.12-slim AS data

# Either 'GRCh38' or 'GRCh37'
ARG ASSEMBLY="GRCh38"

# Tell build-seqrepo where to put the data
ENV SEQREPO_ROOT_DIR=/seqrepo-${ASSEMBLY}

# Install curl for downloading
RUN apt-get update && apt-get install -y curl

WORKDIR /data

# Download the appropriate genomic reference file based on assembly
RUN if [ "$ASSEMBLY" = "GRCh38" ]; then \
        curl -O https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz; \
    elif [ "$ASSEMBLY" = "GRCh37" ]; then \
        curl -O https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.fna.gz; \
    else \
        echo "Unknown assembly: $ASSEMBLY" && exit 1; \
    fi

# Builder image
FROM python:3.12-slim AS build

# Either 'GRCh38' or 'GRCh37'
ARG ASSEMBLY="GRCh38"

# Install packages needed for the build
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    curl \
    git \
    libpq-dev \
    python3-pip \
    python3-venv \
    tabix \
    rsync \
    zlib1g-dev \
    postgresql \
    unzip \
    libhts3 \
    ;

WORKDIR /vrs-python

# Copy downloaded genomic files from data layer
COPY --from=data /data/*.fna.gz /vrs-python/

# Setup the virtual env for vrs-python
RUN python3 -m venv /vrs-python/venv
ENV PATH=/vrs-python/venv/bin:$PATH

# Tell build-seqrepo where to put the data
ENV SEQREPO_ROOT_DIR=/seqrepo-${ASSEMBLY}

# Install vrs-python
RUN /vrs-python/venv/bin/python3 -m pip install -U setuptools
RUN /vrs-python/venv/bin/python3 -m pip install 'ga4gh.vrs[extras]' biocommons.seqrepo

COPY build-${ASSEMBLY}.bash /vrs-python/build-seqrepo.bash
RUN bash /vrs-python/build-seqrepo.bash

# Final image
FROM python:3.12-slim AS vrs-python
ARG ASSEMBLY
ENV ASSEMBLY=${ASSEMBLY}

# Install runtime required packages
RUN apt-get update && apt-get install -y libpq-dev

# Copy over artifacts from the builder
COPY --from=build /vrs-python /vrs-python
COPY --from=build /seqrepo-${ASSEMBLY} /seqrepo-${ASSEMBLY}

# Copy over run script
COPY ./entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Set environment variables
ENV GA4GH_VRS_DATAPROXY_URI="seqrepo+file:///seqrepo-${ASSEMBLY}/master"
ENV SEQREPO_ROOT_DIR=/seqrepo-${ASSEMBLY}
ENV VIRTUAL_ENV=/vrs-python/venv
ENV PATH=/vrs-python/venv/bin:$PATH

WORKDIR /

ENTRYPOINT [ "/entrypoint.sh" ]
