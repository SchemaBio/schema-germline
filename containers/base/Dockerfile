# =============================================================================
# Schema Germline - Base Environment
# =============================================================================
# Purpose: Minimal Python environment for auxiliary scripts and glue code
# Note: Bioinformatics tools (samtools, bcftools, etc.) are NOT included here
#       They are handled by specialized containers to avoid duplication
# =============================================================================

# -----------------------------------------------------------------------------
# Base Image
# -----------------------------------------------------------------------------
FROM ubuntu:22.04 AS base

LABEL maintainer="schemabio"
LABEL description="Base environment for schema-germline pipeline glue scripts"
LABEL version="2.0.0"

# -----------------------------------------------------------------------------
# Environment
# -----------------------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONUNBUFFERED=1 \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# -----------------------------------------------------------------------------
# Install Minimal Dependencies
# -----------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Core utilities
    wget \
    curl \
    git \
    tzdata \
    ca-certificates \
    # File processing
    findutils \
    coreutils \
    sed \
    awk \
    grep \
    # Compression
    gzip \
    bzip2 \
    xz-utils \
    zip \
    unzip \
    # Python
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    # Data formats
    jq \
    yq \
    # Other utilities
    procps \
    vim-tiny \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# -----------------------------------------------------------------------------
# Python Virtual Environment
# -----------------------------------------------------------------------------
RUN python3 -m venv /opt/venv && \
    . /opt/venv/bin/activate && \
    pip install --upgrade pip wheel setuptools

# -----------------------------------------------------------------------------
# Install Python Packages
# -----------------------------------------------------------------------------
COPY containers/base/requirements.txt /tmp/requirements.txt

RUN . /opt/venv/bin/activate && \
    pip install -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

# -----------------------------------------------------------------------------
# Copy Glue Scripts
# -----------------------------------------------------------------------------
COPY bin/ /opt/scripts/
RUN chmod +x /opt/scripts/*.py /opt/scripts/*.sh 2>/dev/null || true

ENV PATH="/opt/scripts:${PATH}"

# -----------------------------------------------------------------------------
# Working Directory
# -----------------------------------------------------------------------------
WORKDIR /workspace

# -----------------------------------------------------------------------------
# Entrypoint
# -----------------------------------------------------------------------------
ENTRYPOINT ["/bin/bash", "-c"]
