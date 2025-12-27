# Schema Germline Pipeline

基于 Nextflow DSL2 的胚系变异检测流程。

## 快速开始

### 1. 环境要求

- Nextflow >= 23.04.0
- Docker / Singularity / Podman (任选其一)

### 2. 准备输入文件

**样本表 (samplesheet.csv)**
```csv
sample_id,read1,read2
sample1,/path/to/sample1_R1.fq.gz,/path/to/sample1_R2.fq.gz
sample2,/path/to/sample2_R1.fq.gz,/path/to/sample2_R2.fq.gz
```

**参考基因组**
- FASTA 文件 (需要 `.fai` 索引)
- BWA-MEM2 或 BWA 索引目录

### 3. 运行流程

```bash
# 本地 Docker 运行
nextflow run main.nf \
    --input samplesheet.csv \
    --fasta /path/to/reference.fa \
    --bwamem2_index /path/to/bwamem2_index \
    --outdir ./results \
    -profile docker

# HPC Singularity 运行
nextflow run main.nf \
    --input samplesheet.csv \
    --fasta /path/to/reference.fa \
    --bwamem2_index /path/to/bwamem2_index \
    --outdir ./results \
    -profile slurm,singularity

# 低内存环境 (使用 bwa 替代 bwa-mem2)
nextflow run main.nf \
    --input samplesheet.csv \
    --fasta /path/to/reference.fa \
    --bwa_index /path/to/bwa_index \
    --outdir ./results \
    -profile docker
```

### 4. 主要参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--input` | 样本表 CSV 文件 | 必填 |
| `--fasta` | 参考基因组 FASTA | 必填 |
| `--fasta_fai` | FASTA 索引文件 | `${fasta}.fai` |
| `--bwamem2_index` | BWA-MEM2 索引目录 | - |
| `--bwa_index` | BWA 索引目录 (低内存备选) | - |
| `--outdir` | 输出目录 | `./results` |
| `--save_trimmed_reads` | 保存过滤后的 FASTQ | `false` |
| `--save_markdup_bam` | 保存去重后的 BAM | `false` |

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
