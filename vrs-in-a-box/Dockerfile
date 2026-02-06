# Builder image
FROM python:3.12-slim AS build

# Either 'GRCh38' or 'GRCh37'
ARG ASSEMBLY="GRCh38"
ARG VRS_PYTHON_VERSION="2.1.2"

# Install packages needed for the build
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    curl \
    git \
    libpq-dev \
    python3-pip \
    python3-venv \
    rsync \
    zlib1g-dev \
    postgresql \
    unzip \
    ;

WORKDIR /vrs-python

# Setup the virtual env for vrs-python
RUN python3 -m venv /vrs-python/venv
ENV PATH=/vrs-python/venv/bin:$PATH

# Install vrs-python
RUN /vrs-python/venv/bin/python3 -m pip install -U setuptools
RUN /vrs-python/venv/bin/python3 -m pip install ga4gh.vrs[extras]==${VRS_PYTHON_VERSION}

# Download and unpack seqrepo files
COPY /seqrepo-${ASSEMBLY}.zip /
RUN unzip /seqrepo-${ASSEMBLY}.zip -d /

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
COPY ./run.sh /run.sh

# Set environment variables
ENV GA4GH_VRS_DATAPROXY_URI="seqrepo+file:///seqrepo-${ASSEMBLY}/master"
ENV VIRTUAL_ENV=/vrs-python/venv
ENV PATH=/vrs-python/venv/bin:$PATH

WORKDIR /

ENTRYPOINT [ "/run.sh" ]
