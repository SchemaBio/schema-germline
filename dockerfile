FROM ubuntu:22.04

LABEL maintainer="SchemaBio"

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /opt

# 安装 Python 和你要求的工具
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    bcftools \
    bedtools \
    tabix \
    git \
    && rm -rf /var/lib/apt/lists/*

# 安装 Python 库
RUN pip3 install --break-system-packages pysam pandas openpyxl polars

# 拉取项目
RUN git clone https://github.com/schemabio/schema-germline.git /opt/schema-germline

WORKDIR /opt/schema-germline

CMD ["/bin/bash"]