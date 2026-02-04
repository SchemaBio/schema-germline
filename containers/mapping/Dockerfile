# =============================================================================
# Schema Germline - Mapping Container
# =============================================================================
# Purpose: Sequence alignment and QC tools
# Tools: samtools, bcftools, bwa, bwa-mem2, fastp, mosdepth, bamdst, xamdst
# Note: This image contains core bioinformatics tools used in multiple steps
# =============================================================================

# -----------------------------------------------------------------------------
# Base Image
# -----------------------------------------------------------------------------
FROM ubuntu:22.04

LABEL maintainer="schemabio"
LABEL description="Mapping tools for schema-germline pipeline"
LABEL version="2.0.0"

# -----------------------------------------------------------------------------
# Version Definitions
# -----------------------------------------------------------------------------
ARG SAMTOOLS_VERSION=1.21
ARG BCFTOOLS_VERSION=1.21
ARG BWA_VERSION=0.7.19
ARG BWA_MEM2_VERSION=2.2.1
ARG FASTP_VERSION=0.24.2
ARG MOSDEPTH_VERSION=0.7.1
ARG Sambamba_VERSION=1.0.2

ENV DEBIAN_FRONTEND=noninteractive

# -----------------------------------------------------------------------------
# Install Build Dependencies
# -----------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    make \
    cmake \
    autoconf \
    automake \
    libtool \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libdeflate-dev \
    libncurses-dev \
    wget \
    curl \
    git \
    ca-certificates \
    unzip \
    pkg-config \
    libhts-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /software

# =============================================================================
# samtools (核心工具 - 会在多个镜像中重复)
# =============================================================================
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local --with-htslib=system && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf samtools-${SAMTOOLS_VERSION}*

# =============================================================================
# bcftools (核心工具 - 会在多个镜像中重复)
# =============================================================================
RUN wget -q https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --prefix=/usr/local --with-htslib=system && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf bcftools-${BCFTOOLS_VERSION}*

# =============================================================================
# htslib (samtools/bcftools 依赖)
# =============================================================================
RUN wget -q https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd htslib-${SAMTOOLS_VERSION} && \
    autoreconf -i && \
    ./configure --prefix=/usr/local --with-htslib=system && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf htslib-${SAMTOOLS_VERSION}*

# =============================================================================
# bwa
# =============================================================================
RUN wget -q https://github.com/lh3/bwa/archive/refs/tags/v${BWA_VERSION}.tar.gz && \
    tar -xzf v${BWA_VERSION}.tar.gz && \
    cd bwa-${BWA_VERSION} && \
    make -j$(nproc) && \
    cp bwa /usr/local/bin/ && \
    cd .. && \
    rm -rf bwa-${BWA_VERSION} v${BWA_VERSION}.tar.gz

# =============================================================================
# bwa-mem2
# =============================================================================
RUN wget -q https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_MEM2_VERSION}/bwa-mem2-${BWA_MEM2_VERSION}_x64-linux.tar.bz2 && \
    tar -xjf bwa-mem2-${BWA_MEM2_VERSION}_x64-linux.tar.bz2 && \
    cp bwa-mem2-${BWA_MEM2_VERSION}_x64-linux/bwa-mem2* /usr/local/bin/ && \
    chmod 755 /usr/local/bin/bwa-mem2* && \
    rm -rf bwa-mem2-${BWA_MEM2_VERSION}_x64-linux bwa-mem2-${BWA_MEM2_VERSION}_x64-linux.tar.bz2

# =============================================================================
# fastp
# =============================================================================
RUN wget -q http://opengene.org/fastp/fastp && \
    chmod a+x fastp && \
    mv fastp /usr/local/bin/

# =============================================================================
# mosdepth
# =============================================================================
RUN wget -q https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth && \
    chmod a+x mosdepth && \
    mv mosdepth /usr/local/bin/

# =============================================================================
# mosdepth d4 (optional, for coverage in specific formats)
# =============================================================================
RUN wget -q https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth_d4 && \
    chmod a+x mosdepth_d4 && \
    mv mosdepth_d4 /usr/local/bin/

# =============================================================================
# sambamba (备份去重工具)
# =============================================================================
RUN wget -q https://github.com/biod/sambamba/releases/download/v${Sambamba_VERSION}/sambamba-${Sambamba_VERSION}-linux-amd64-static.gz && \
    gunzip sambamba-${Sambamba_VERSION}-linux-amd64-static.gz && \
    mv sambamba-${Sambamba_VERSION}-linux-amd64-static /usr/local/bin/sambamba && \
    chmod 755 /usr/local/bin/sambamba

# =============================================================================
# bamdst (覆盖度统计)
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends libhts-dev && \
    rm -rf /var/lib/apt/lists/* && \
    git clone --depth 1 https://github.com/shiquan/bamdst.git && \
    cd bamdst && \
    make && \
    cp bamdst /usr/local/bin/ && \
    cd .. && \
    rm -rf bamdst

# =============================================================================
# xamdst (增强覆盖度统计)
# =============================================================================
RUN git clone --depth 1 https://github.com/pzweuj/xamdst.git && \
    cd xamdst && \
    make && \
    cp xamdst /usr/local/bin/ && \
    cd .. && \
    rm -rf xamdst

# =============================================================================
# tabix + bgzip (VCF 索引工具)
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    tabix \
    && rm -rf /var/lib/apt/lists/*

# =============================================================================
# Python + pybedtools (用于脚本)
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir pybedtools

# =============================================================================
# Cleanup
# =============================================================================
RUN apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /software/*

# =============================================================================
# Environment
# =============================================================================
ENV PATH="/usr/local/bin:${PATH}"
WORKDIR /data
CMD ["/bin/bash"]
