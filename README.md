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
  "read1": "/data/sample1_R1.fq.gz",
  "read2": "/data/sample1_R2.fq.gz",
  "reference": {
    "fasta": "/reference/hg38.fa",
    "fasta_fai": "/reference/hg38.fa.fai",
    "bwamem2_index": "/reference/bwamem2_index"
  },
  "outdir": "./results"
}
```

> 注意：路径应为容器内路径，与 `-v` 挂载目录对应

### 3. 运行流程

```bash
# 运行单个样本
docker run --rm -it \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v $(pwd):/workspace \
    -v /path/to/data:/data \
    -v /path/to/reference:/reference \
    -w /workspace \
    nextflow/nextflow:latest \
    nextflow run main.nf \
        --config sample1.json \
        -profile docker

# 批量运行多个样本
for config in samples/*.json; do
    docker run --rm -d \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v $(pwd):/workspace \
        -v /path/to/data:/data \
        -v /path/to/reference:/reference \
        -w /workspace \
        nextflow/nextflow:latest \
        nextflow run main.nf \
            --config $config \
            -profile docker
done

# HPC Singularity 运行
nextflow run main.nf \
    --config sample1.json \
    -profile slurm,singularity
```

### 4. 配置文件说明

```json
{
  "sample_id": "样本ID",
  "read1": "R1 fastq 路径",
  "read2": "R2 fastq 路径",
  "reference": {
    "fasta": "参考基因组 FASTA",
    "fasta_fai": "FASTA 索引 (可选)",
    "bwamem2_index": "BWA-MEM2 索引目录",
    "bwa_index": "BWA 索引目录 (低内存备选)"
  },
  "outdir": "输出目录 (默认 ./results)"
}
```

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
