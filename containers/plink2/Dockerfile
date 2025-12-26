# 1. 锁定上游官方版本
FROM alpine:3.23.2

# 2. 注入绳墨生物的"基因"
LABEL maintainer="Schema Bio"
LABEL vendor="Schema Bio"
LABEL schema.version="20251226"

# 3. 版本号定义
ARG PLINK2_DATE=20251205

# 4. 安装必要依赖
RUN apk update && \
    apk add --no-cache \
    bash \
    wget \
    unzip \
    libc6-compat \
    libstdc++

WORKDIR /software

# 5. 下载三个版本的 plink2
# Intel AVX2 (现代 Intel/AMD 处理器)
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_${PLINK2_DATE}.zip && \
    unzip plink2_linux_avx2_${PLINK2_DATE}.zip -d plink2_avx2 && \
    mv plink2_avx2/plink2 /usr/local/bin/plink2_avx2 && \
    chmod 755 /usr/local/bin/plink2_avx2 && \
    rm -rf plink2_linux_avx2_${PLINK2_DATE}.zip plink2_avx2

# AMD AVX2 (AMD Zen 架构优化)
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_amd_avx2_${PLINK2_DATE}.zip && \
    unzip plink2_linux_amd_avx2_${PLINK2_DATE}.zip -d plink2_amd && \
    mv plink2_amd/plink2 /usr/local/bin/plink2_amd && \
    chmod 755 /usr/local/bin/plink2_amd && \
    rm -rf plink2_linux_amd_avx2_${PLINK2_DATE}.zip plink2_amd

# x86_64 通用版本 (兼容性最好)
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_${PLINK2_DATE}.zip && \
    unzip plink2_linux_x86_64_${PLINK2_DATE}.zip -d plink2_x86 && \
    mv plink2_x86/plink2 /usr/local/bin/plink2_x86_64 && \
    chmod 755 /usr/local/bin/plink2_x86_64 && \
    rm -rf plink2_linux_x86_64_${PLINK2_DATE}.zip plink2_x86

# 6. 创建自动选择脚本
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

# 7. 清理
RUN rm -rf /software/*

WORKDIR /data
ENV PATH="/usr/local/bin:${PATH}"
CMD ["/bin/bash"]
