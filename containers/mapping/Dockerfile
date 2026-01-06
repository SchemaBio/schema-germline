# 1. 锁定上游官方版本
FROM ubuntu:22.04

# 2. 注入绳墨生物的"基因"
LABEL maintainer="Schema Bio"
LABEL vendor="Schema Bio"
LABEL schema.version="20260106"

# 3. 版本号定义 - 方便后续更新
ARG SAMTOOLS_VERSION=1.23
ARG BCFTOOLS_VERSION=1.23
ARG BWA_VERSION=0.7.19
ARG SAMBAMBA_VERSION=1.0.1
ARG BWA_MEM2_VERSION=2.3
ARG MOSDEPTH_VERSION=0.3.12

# Install necessary packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    bash \
    vim \
    ca-certificates \
    zlib1g-dev \
    libbz2-dev \
    wget \
    unzip \
    build-essential \
    libncurses-dev \
    xz-utils \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libdeflate-dev \
    libtool \
    autoconf \
    automake \
    git

WORKDIR /software

# samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    rm samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix $(pwd) && \
    make && \
    cp samtools /usr/local/bin/ && \
    rm -rf /software/samtools-${SAMTOOLS_VERSION}

# bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -jxvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    rm bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --prefix $(pwd) && \
    make && \
    cp bcftools /usr/local/bin/ && \
    rm -rf /software/bcftools-${BCFTOOLS_VERSION}

# bwa
RUN wget https://github.com/lh3/bwa/archive/refs/tags/v${BWA_VERSION}.zip && \
    unzip v${BWA_VERSION}.zip && \
    cd bwa-${BWA_VERSION} && \
    make && \
    cp bwa /usr/local/bin/ && \
    cd .. && \
    rm -rf /software/v${BWA_VERSION}.zip && \
    rm -rf /software/bwa-${BWA_VERSION}

# sambamba
RUN wget https://github.com/biod/sambamba/releases/download/v${SAMBAMBA_VERSION}/sambamba-${SAMBAMBA_VERSION}-linux-amd64-static.gz && \
    gunzip sambamba-${SAMBAMBA_VERSION}-linux-amd64-static.gz && \
    mv sambamba-${SAMBAMBA_VERSION}-linux-amd64-static /usr/local/bin/ && \
    chmod 777 /usr/local/bin/sambamba-${SAMBAMBA_VERSION}-linux-amd64-static && \
    ln -s /usr/local/bin/sambamba-${SAMBAMBA_VERSION}-linux-amd64-static /usr/local/bin/sambamba

# bamdst
RUN wget https://github.com/shiquan/bamdst/archive/refs/heads/master.zip && \
    unzip master.zip && \
    cd bamdst-master && \
    make && \
    cp bamdst /usr/local/bin/ && \
    rm -rf /software/master.zip && \
    rm -rf /software/bamdst-master

# xamdst
RUN apt-get install -y libhts-dev pkg-config && \
    git clone https://github.com/pzweuj/xamdst.git && \
    cd xamdst && \
    git pull && \
    make && \
    cp xamdst /usr/local/bin/ && \
    rm -rf /software/xamdst

# fastp - prebuilt binary
RUN wget -O /usr/local/bin/fastp http://opengene.org/fastp/fastp && \
    chmod a+x /usr/local/bin/fastp

# bwa-mem2
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_MEM2_VERSION}/bwa-mem2-${BWA_MEM2_VERSION}_x64-linux.tar.bz2 && \
    tar -jxvf bwa-mem2-${BWA_MEM2_VERSION}_x64-linux.tar.bz2 && \
    mv bwa-mem2-${BWA_MEM2_VERSION}_x64-linux/bwa-mem2* /usr/local/bin/ && \
    chmod 777 /usr/local/bin/bwa-mem2*

# mosdepth D4
RUN wget https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth_d4 && \
    mv mosdepth_d4 /usr/local/bin/ && \
    chmod 777 /usr/local/bin/mosdepth_d4

# mosdepth
RUN wget https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth && \
    mv mosdepth /usr/local/bin/ && \
    chmod 777 /usr/local/bin/mosdepth

# default command
WORKDIR /usr/local/bin/
ENV PATH="/usr/local/bin/:${PATH}"
CMD ["/bin/bash"]
