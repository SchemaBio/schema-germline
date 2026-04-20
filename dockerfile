# schema-germline Pipeline Docker Image
# 包含整个流程及 bin 脚本所需的 Python 库

FROM python:3.11-slim-bookworm

LABEL maintainer="SchemaBio"
LABEL description="Germline variant analysis pipeline"

# 设置工作目录
WORKDIR /pipeline

# 安装系统依赖 (生物信息工具)
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    bedtools \
    bcftools \
    tabix \
    && rm -rf /var/lib/apt/lists/*

# 安装 Python 依赖 (bin 脚本所需)
# make_pedigree.py: PyYAML
# results_to_parquet.py: pandas, pyarrow
# cnv_segment.py: numpy, scipy
RUN pip install --no-cache-dir \
    PyYAML \
    pandas \
    pyarrow \
    numpy \
    scipy

# 复制项目文件到 /pipeline (排除 .git)
COPY bin/ /pipeline/bin/
COPY conf/ /pipeline/conf/
COPY modules/ /pipeline/modules/
COPY workflows/ /pipeline/workflows/
COPY assets/ /pipeline/assets/
COPY examples/ /pipeline/examples/
COPY docs/ /pipeline/docs/
COPY main.nf /pipeline/main.nf
COPY nextflow.config /pipeline/nextflow.config
COPY README.md /pipeline/README.md
COPY LICENSE /pipeline/LICENSE

# 设置 PATH 包含 bin 目录
ENV PATH="/pipeline/bin:${PATH}"

# 设置默认入口
ENTRYPOINT ["/bin/bash"]