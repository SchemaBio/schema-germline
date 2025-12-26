# 第一阶段：构建
FROM python:3.9-slim AS builder

# 注入绳墨生物的"基因"
LABEL maintainer="Schema Bio"
LABEL vendor="Schema Bio"
LABEL schema.version="20251226"

# 安装必要的依赖
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        gcc zlib1g-dev libhts-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# 安装 Python 依赖
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir whatshap

# 第二阶段：最终镜像
FROM python:3.9-slim

# 从第一阶段复制必要的文件
COPY --from=builder /usr/local/lib/python3.9/site-packages /usr/local/lib/python3.9/site-packages
COPY --from=builder /usr/local/bin/whatshap /usr/local/bin/whatshap

# 安装运行时依赖
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        zlib1g libhts3 tabix && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# 默认命令
CMD ["whatshap", "--help"]