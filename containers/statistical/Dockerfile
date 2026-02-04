# =============================================================================
# Schema Germline - Statistical Container
# =============================================================================
# Purpose: Statistical analysis and population genetics tools
# Tools: Plink2, Peddy, Vcftools, Beagle (optional)
# Note: samtools is included for VCF/BCF operations
# =============================================================================

# -----------------------------------------------------------------------------
# Base Image
# -----------------------------------------------------------------------------
FROM alpine:3.20

LABEL maintainer="schemabio"
LABEL description="Statistical analysis tools for schema-germline pipeline"
LABEL version="2.0.0"

# -----------------------------------------------------------------------------
# Version Definitions
# -----------------------------------------------------------------------------
ARG PLINK2_DATE=20251205
ARG PEDDY_VERSION=0.4.3
ARG SAMTOOLS_VERSION=1.21

ENV PLINK2_VERSION=${PLINK2_DATE}

# =============================================================================
# System Dependencies
# =============================================================================
RUN apk add --no-cache \
    # Build tools
    build-base \
    gcc \
    g++ \
    make \
    cmake \
    autoconf \
    automake \
    libtool \
    # samtools dependencies
    zlib-dev \
    bzip2-dev \
    xz-dev \
    curl-dev \
    openssl-dev \
    # Python
    python3 \
    py3-pip \
    py3-numpy \
    # Other tools
    wget \
    bash \
    coreutils \
    unzip \
    tar \
    libstdc++ \
    libc6-compat

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
# Plink2 (多版本支持 CPU 架构优化)
# =============================================================================
# Intel AVX2 (现代 Intel/AMD 处理器)
RUN wget -q https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_${PLINK2_DATE}.zip && \
    unzip -q plink2_linux_avx2_${PLINK2_DATE}.zip -d plink2_avx2 && \
    mv plink2_avx2/plink2 /usr/local/bin/plink2_avx2 && \
    chmod 755 /usr/local/bin/plink2_avx2 && \
    rm -rf plink2_linux_avx2_${PLINK2_DATE}.zip plink2_avx2

# AMD AVX2 (AMD Zen 架构优化)
RUN wget -q https://s3.amazonaws.com/plink2-assets/plink2_linux_amd_avx2_${PLINK2_DATE}.zip && \
    unzip -q plink2_linux_amd_avx2_${PLINK2_DATE}.zip -d plink2_amd && \
    mv plink2_amd/plink2 /usr/local/bin/plink2_amd && \
    chmod 755 /usr/local/bin/plink2_amd && \
    rm -rf plink2_linux_amd_avx2_${PLINK2_DATE}.zip plink2_amd

# x86_64 通用版本 (兼容性最好)
RUN wget -q https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_${PLINK2_DATE}.zip && \
    unzip -q plink2_linux_x86_64_${PLINK2_DATE}.zip -d plink2_x86 && \
    mv plink2_x86/plink2 /usr/local/bin/plink2_x86_64 && \
    chmod 755 /usr/local/bin/plink2_x86_64 && \
    rm -rf plink2_linux_x86_64_${PLINK2_DATE}.zip plink2_x86

# 自动选择脚本
RUN echo '#!/bin/bash' > /usr/local/bin/plink2 && \
    echo '# 自动检测 CPU 并选择最优版本' >> /usr/local/bin/plink2 && \
    echo 'if grep -q "AuthenticAMD" /proc/cpuinfo && grep -q "avx2" /proc/cpuinfo; then' >> /usr/local/bin/plink2 && \
    echo '    exec /usr/local/bin/plink2_amd "$@"' >> /usr/local/bin/plink2 && \
    echo 'elif grep -q "avx2" /proc/cpuinfo; then' >> /usr/local/bin/plink2 && \
    echo '    exec /usr/local/bin/plink2_avx2 "$@"' >> /usr/local/bin/plink2 && \
    echo 'else' >> /usr/local/bin/plink2 && \
    echo '    exec /usr/local/bin/plink2_x86_64 "$@"' >> /usr/local/bin/plink2 && \
    echo 'fi' >> /usr/local/bin/plink2 && \
    chmod 755 /usr/local/bin/plink2

# =============================================================================
# Peddy (家系和群体结构检查)
# =============================================================================
RUN pip3 install --no-cache-dir peddy==${PEDDY_VERSION}

# =============================================================================
# Vcftools (VCF 文件操作)
# =============================================================================
RUN apk add --no-cache vcftools || true

# =============================================================================
# Beagle (基因型填充 - 可选)
# =============================================================================
ARG BEAGLE_VERSION=5.4
RUN wget -q https://faculty.washington.edu/browning/beagle/${BEAGLE_VERSION}/beagle.jar -O /usr/local/bin/beagle.jar && \
    chmod 644 /usr/local/bin/beagle.jar

# =============================================================================
# Python packages for statistics
# =============================================================================
RUN pip3 install --no-cache-dir \
    pandas \
    scipy \
    matplotlib \
    seaborn

# =============================================================================
# Cleanup
# =============================================================================
RUN rm -rf /software/* /var/cache/apk/* /tmp/*

# =============================================================================
# Environment
# =============================================================================
ENV PATH="/usr/local/bin:${PATH}" \
    LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

WORKDIR /data
CMD ["plink2"]
