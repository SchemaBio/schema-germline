# =============================================================================
# Schema Germline - CNV Container
# =============================================================================
# Purpose: Copy Number Variation analysis tools
# Tools: CNVkit, Stranger, SVDB, bedtools
# Note: samtools is included for VCF/BCF operations
# =============================================================================

# -----------------------------------------------------------------------------
# Base Image
# -----------------------------------------------------------------------------
FROM python:3.11-slim

LABEL maintainer="schemabio"
LABEL description="CNV analysis tools for schema-germline pipeline"
LABEL version="2.0.0"

# -----------------------------------------------------------------------------
# Version Definitions
# -----------------------------------------------------------------------------
ARG CNVKIT_VERSION=0.9.10
ARG STRANGER_VERSION=0.8.0
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
    wget \
    curl \
    git \
    ca-certificates \
    # samtools dependencies
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libdeflate-dev \
    libncurses-dev \
    # Other tools
    tabix \
    unzip \
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
RUN apt-get update && apt-get install -y --no-install-recommends \
    autoconf \
    automake \
    libtool \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q https://github.com/samtools/bcftools/releases/download/${SAMTOOLS_VERSION}/bcftools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf bcftools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local --with-htslib=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf bcftools-${SAMTOOLS_VERSION}*

# =============================================================================
# Python Dependencies for CNVkit
# =============================================================================
RUN pip install --no-cache-dir --upgrade pip wheel setuptools && \
    pip install --no-cache-dir \
    numpy \
    scipy \
    pandas \
    pysam \
    pybedtools

# =============================================================================
# CNVkit
# =============================================================================
RUN pip install --no-cache-dir cnvkit==${CNVKIT_VERSION}

# =============================================================================
# Stranger (STR 重复次数注释)
# =============================================================================
RUN pip install --no-cache-dir stranger==${STRANGER_VERSION}

# =============================================================================
# SVDB (结构变异数据库注释和合并)
# =============================================================================
RUN pip install --no-cache-dir svdb

# =============================================================================
# bedtools (基因组区间操作)
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# =============================================================================
# R + DNAcopy (可选：用于分割算法)
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('DNAcopy', repos='https://cloud.r-project.org', quiet=TRUE)"

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
CMD ["cnvkit"]
