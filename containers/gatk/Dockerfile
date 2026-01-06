FROM broadinstitute/gatk:4.6.2.0

# 注入绳墨生物的"基因"
LABEL maintainer="Schema Bio"
LABEL vendor="Schema Bio"
LABEL schema.version="20251226"

# 安装 samtools 用于 CRAM 索引
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    && rm -rf /var/lib/apt/lists/*