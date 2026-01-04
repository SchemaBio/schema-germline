# =============================================================================
# Schema Germline - Base Environment for Glue Scripts
# =============================================================================
# Purpose: Minimal environment for running auxiliary Python/R scripts
# Note: Bioinformatics tools (samtools, bcftools, etc.) are handled by
#       Nextflow module containers, not this base image
# =============================================================================

# -----------------------------------------------------------------------------
# Base Image
# -----------------------------------------------------------------------------
FROM ubuntu:22.04 AS base

LABEL maintainer="Schema Bio"
LABEL description="Base environment for schema-germline pipeline glue scripts"
LABEL version="1.0.0"

# -----------------------------------------------------------------------------
# Environment
# -----------------------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONUNBUFFERED=1

# -----------------------------------------------------------------------------
# Install Minimal Dependencies
# -----------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Core utilities for scripts
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
    # Data formats (for config parsing)
    jq \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# -----------------------------------------------------------------------------
# Python Virtual Environment
# -----------------------------------------------------------------------------
RUN python3 -m venv /opt/venv && \
    . /opt/venv/bin/activate && \
    pip install --upgrade pip

# -----------------------------------------------------------------------------
# Install Python Packages
# ----------------------------------------------------------------------------=
COPY requirements.txt /tmp/requirements.txt

RUN . /opt/venv/bin/activate && \
    pip install -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

# -----------------------------------------------------------------------------
# Copy Glue Scripts from project bin/
# -----------------------------------------------------------------------------
# Context: Build from project root, e.g.:
#   docker build -f containers/base/Dockerfile -t base .
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
