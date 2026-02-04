# =============================================================================
# Schema Germline - SV Container
# =============================================================================
# Purpose: Structural Variation detection tools
# Tools: ExpansionHunter, Surf, Manta (optional)
# Note: samtools is included for BAM/CRAM operations
# =============================================================================

# -----------------------------------------------------------------------------
# Base Image
# -----------------------------------------------------------------------------
FROM ubuntu:22.04

LABEL maintainer="schemabio"
LABEL description="SV detection tools for schema-germline pipeline"
LABEL version="2.0.0"

# -----------------------------------------------------------------------------
# Version Definitions
# -----------------------------------------------------------------------------
ARG EXPANSIONHUNTER_VERSION=5.0.0
ARG SURF_VERSION=2.0
ARG SAMTOOLS_VERSION=1.21

ENV DEBIAN_FRONTEND=noninteractive

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
    libboost-all-dev \
    # Other tools
    wget \
    curl \
    git \
    ca-certificates \
    unzip \
    tabix \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /software

# =============================================================================
# samtools (用于 BAM/CRAM 操作)
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
# ExpansionHunter (STR 检测)
# =============================================================================
RUN wget -q https://github.com/Illumina/ExpansionHunterDenovo/releases/download/v${EXPANSIONHUNTER_VERSION}/ExpansionHunterDenovo-v${EXPANSIONHUNTER_VERSION}-linux.tar.gz && \
    tar -xzf ExpansionHunterDenovo-v${EXPANSIONHUNTER_VERSION}-linux.tar.gz && \
    cp ExpansionHunterDenovo-v${EXPANSIONHUNTER_VERSION}-linux/ExpansionHunter* /usr/local/bin/ && \
    chmod 755 /usr/local/bin/ExpansionHunter* && \
    rm -rf ExpansionHunterDenovo-v${EXPANSIONHUNTER_VERSION}-linux*

# =============================================================================
# Surf2 (STR 检测增强)
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir surf2==${SURF_VERSION}

# =============================================================================
# Manta (结构变异检测 - 可选，用于全基因组 SV)
# =============================================================================
ARG MANTA_VERSION=1.6.0
RUN wget -q https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.centos6_x86_64.tar.gz && \
    tar -xzf manta-${MANTA_VERSION}.centos6_x86_64.tar.gz && \
    cp manta-${MANTA_VERSION}.centos6_x86_64/bin/manta* /usr/local/bin/ && \
    cp -r manta-${MANTA_VERSION}.centos6_x86_64/lib /usr/local/ && \
    chmod 755 /usr/local/bin/manta* && \
    rm -rf manta-${MANTA_VERSION}.centos6_x86_64*

# =============================================================================
# Python + Biopython (用于脚本)
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir biopython pybedtools

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
CMD ["/bin/bash"]
