# =============================================================================
# Schema Germline - Core Container
# =============================================================================
# Purpose: Core variant analysis and filtering tools
# Tools: GenMod, Slivar, WhatsHap, BCFtools
# Note: samtools is included for VCF/BCF operations
# =============================================================================

# -----------------------------------------------------------------------------
# Base Image
# -----------------------------------------------------------------------------
FROM python:3.11-slim

LABEL maintainer="schemabio"
LABEL description="Core analysis tools for schema-germline pipeline"
LABEL version="2.0.0"

# -----------------------------------------------------------------------------
# Version Definitions
# -----------------------------------------------------------------------------
ARG GENMOD_VERSION=3.10.2
ARG SLIVAR_VERSION=3.0.5
ARG WHATSHAP_VERSION=2.3
ARG SAMTOOLS_VERSION=1.21

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1

# =============================================================================
# System Dependencies
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Build tools
    build-essential \
    gcc \
    g++ \
    make \
    cmake \
    autoconf \
    automake \
    libtool \
    # samtools dependencies
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libdeflate-dev \
    libncurses-dev \
    # Python dependencies
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    # Other tools
    wget \
    curl \
    git \
    ca-certificates \
    unzip \
    tabix \
    libffi-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /software

# =============================================================================
# samtools (用于 VCF/BCF 操作)
# =============================================================================
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf samtools-${SAMTOOLS_VERSION}*

# =============================================================================
# htslib (samtools 依赖)
# =============================================================================
RUN wget -q https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd htslib-${SAMTOOLS_VERSION} && \
    autoreconf -i && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf htslib-${SAMTOOLS_VERSION}*

# =============================================================================
# bcftools (用于 VCF 操作)
# =============================================================================
RUN wget -q https://github.com/samtools/bcftools/releases/download/${SAMTOOLS_VERSION}/bcftools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf bcftools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local --with-htslib=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf bcftools-${SAMTOOLS_VERSION}*

# =============================================================================
# Python Dependencies
# =============================================================================
RUN pip install --no-cache-dir --upgrade pip wheel setuptools

# =============================================================================
# GenMod (遗传模式注释)
# =============================================================================
RUN pip install --no-cache-dir genmod==${GENMOD_VERSION}

# =============================================================================
# Slivar (家系变异过滤)
# =============================================================================
RUN pip install --no-cache-dir slivar==${SLIVAR_VERSION}

# =============================================================================
# WhatsHap (单倍型定相)
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhts-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir whatshap==${WHATSHAP_VERSION}

# =============================================================================
# FreeBayes (可选：用于验证)
# =============================================================================
ARG FREEBAYES_VERSION=1.3.6
RUN apt-get update && apt-get install -y --no-install-recommends \
    tabix \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q https://github.com/freebayes/freebayes/releases/download/v${FREEBAYES_VERSION}/freebayes-${FREEBAYES_VERSION}.tar.gz && \
    tar -xzf freebayes-${FREEBAYES_VERSION}.tar.gz && \
    cd freebayes-${FREEBAYES_VERSION} && \
    make -j$(nproc) && \
    cp bin/freebayes /usr/local/bin/ && \
    cd .. && \
    rm -rf freebayes-${FREEBAYES_VERSION}*

# =============================================================================
# GLNexus (联合分型 - 简化版本)
# =============================================================================
# Note: For production, use the official GLNexus Docker image
# This is a lightweight alternative for basic operations

# =============================================================================
# Additional Python packages
# =============================================================================
RUN pip install --no-cache-dir \
    cyvcf2 \
    numpy \
    pandas \
    scipy

# =============================================================================
# Cleanup
# =============================================================================
RUN apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /software/*

# =============================================================================
# Environment
# =============================================================================
ENV PATH="/usr/local/bin:${PATH}" \
    LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

WORKDIR /data
CMD ["genmod"]
