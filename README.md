# Schema Germline Pipeline

基于 Nextflow DSL2 的胚系变异检测流程。

## 快速开始

### 1. 环境要求

- Docker

### 2. 准备配置文件

每个样本一个 JSON 文件，如 `sample1.json`:

```json
{
  "sample_id": "sample1",
  "read1": "/mnt/d/data/sample1_R1.fq.gz",
  "read2": "/mnt/d/data/sample1_R2.fq.gz",
  "reference": {
    "fasta": "/mnt/d/reference/hg38.fa"
  },
  "outdir": "/mnt/d/analysis/results"
}
```

> bwa/bwa-mem2 索引文件需与 fasta 同目录同前缀，无需单独配置。

### 3. 运行流程

```bash
# 设置路径 (使用真实路径，不要用别名)
WORKDIR=/mnt/d/analysis
DATADIR=/mnt/d/data
REFDIR=/mnt/d/reference

# 运行单个样本
docker run --rm -it \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /path/to/schema-germline:/pipeline:ro \
    -v ${WORKDIR}:${WORKDIR} \
    -v ${DATADIR}:${DATADIR} \
    -v ${REFDIR}:${REFDIR} \
    -w ${WORKDIR} \
    nextflow/nextflow:latest \
    nextflow run /pipeline/main.nf \
        --config sample1.json \
        -profile docker
```

> **Docker-in-Docker 注意事项：** 除 pipeline 目录外，所有路径挂载时容器内外必须一致。因为 Nextflow 调度的子容器由宿主机 Docker daemon 启动，只能看到宿主机文件系统。

### 4. 配置文件说明

```json
{
  "sample_id": "样本ID",
  "read1": "R1 fastq 路径",
  "read2": "R2 fastq 路径",
  "reference": {
    "fasta": "参考基因组 FASTA",
    "fasta_fai": "FASTA 索引 (可选)"
  },
  "outdir": "输出目录 (默认 ./results)"
}
```

> bwa/bwa-mem2 索引文件需与 fasta 同目录同前缀，无需单独配置。

### 5. 输出结构

```
results/
├── fastp/                  # 质控报告 (JSON/HTML)
├── alignment/              # 比对结果 (CRAM/CRAI)
├── markdup/                # 去重统计 (metrics.txt)
└── pipeline_info/          # 流程报告 (timeline/report/trace)
```

## 流程说明

```
FASTP (质控过滤) -> BWA_MEM2 (比对) -> GATK_MARKDUPLICATES (标记重复) -> CRAM
```

## Profiles

| Profile | 说明 |
|---------|------|
| `docker` | 本地 Docker 运行 |
| `singularity` | Singularity 容器 |
| `podman` | Podman 容器 |
| `slurm` | SLURM 调度器 |
| `lsf` | LSF 调度器 |
| `k8s` | Kubernetes |
| `aliyun` | 阿里云 |
| `tencent` | 腾讯云 |

## 项目结构

```
schema-germline/
├── main.nf                 # 入口文件
├── nextflow.config         # 全局配置
├── conf/
│   ├── modules.config      # 模块容器和发布配置
│   ├── aliyun.config       # 阿里云配置
│   └── tencent.config      # 腾讯云配置
├── modules/local/          # 本地模块
│   ├── fastp/main.nf       # 质控过滤
│   ├── bwa_mem2/main.nf    # 比对
│   └── gatk/main.nf        # GATK 工具集
├── workflows/
│   └── wes_single.nf       # WES 单样本流程
└── containers/             # Dockerfile
```

## License

Apache License 2.0
