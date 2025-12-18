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
RUN apt-get update && apt-get install -y curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /data

# Download the appropriate genomic reference file based on assembly
COPY build-${ASSEMBLY}.bash /data/
RUN . /data/build-${ASSEMBLY}.bash \
    && download_reference

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
    postgresql-client \
    unzip \
    libhts3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /vrs-python

# Copy downloaded genomic files from data layer
COPY --from=data /data/* /data/

# Setup the virtual env for vrs-python
RUN python3 -m venv /vrs-python/venv
ENV PATH=/vrs-python/venv/bin:$PATH

# Tell build-seqrepo where to put the data
ENV SEQREPO_ROOT_DIR=/seqrepo-${ASSEMBLY}

# Install vrs-python
RUN /vrs-python/venv/bin/python3 -m pip install -U setuptools 'ga4gh.vrs[extras]' biocommons.seqrepo

# Build the seqrepo data using provided function
COPY build-${ASSEMBLY}.bash /tmp/build-${ASSEMBLY}.bash
RUN cd /data && . /tmp/build-${ASSEMBLY}.bash \
    && build_seqrepo

# Final image
FROM python:3.12-slim AS vrs-python
ARG ASSEMBLY="GRCh38"
ENV ASSEMBLY=${ASSEMBLY}

# Install runtime required packages
RUN apt-get update && apt-get install -y libpq-dev \
    && rm -rf /var/lib/apt/lists/*

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
