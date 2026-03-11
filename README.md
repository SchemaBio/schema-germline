# Schema Germline Pipeline

基于 Nextflow DSL2 的胚系变异检测流程。

## 快速开始

### 1. 环境要求

- Docker
- 自部署镜像已拉取（见下文）

### 2. 拉取镜像

```bash
# Nextflow 镜像
docker pull docker.biotools.space/schemabio/nextflow:25.10.4

# 流程依赖的其他镜像（自动按需拉取）
docker pull docker.biotools.space/schemabio/mapping:2026Jan
docker pull docker.biotools.space/schemabio/gatk:4.6.2.0
docker pull docker.biotools.space/schemabio/cnvkit:v0.9.11.p4
docker pull docker.biotools.space/schemabio/expansionhunter:5.0.0
docker pull docker.biotools.space/schemabio/deepvariant:1.10.0
docker pull docker.biotools.space/schemabio/ensembl-vep:release_115.2
docker pull docker.biotools.space/schemabio/glnexus:v1.4.1
```

### 3. 准备配置文件

每个样本一个 JSON 文件，如 `examples/sample.json`:

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

### 4. 运行流程

#### 完整流程 WES_SINGLE（默认）

```bash
# 设置路径 (使用真实路径，不要用别名)
WORKDIR=/mnt/d/analysis
PIPELINE_DIR=/path/to/schema-germline

docker run --rm -it \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v ${PIPELINE_DIR}:/pipeline:ro \
    -v ${WORKDIR}:${WORKDIR} \
    -w ${WORKDIR} \
    docker.biotools.space/schemabio/nextflow:25.10.4 \
    nextflow run /pipeline/main.nf \
        -config /pipeline/conf/wes_single.config \
        --config /pipeline/examples/sample.json \
        -profile docker
```

#### 质控+比对流程 QC_ALIGNMENT

仅进行 FASTP 质控和 BWA 比对，输出 sorted.cram：

```bash
docker run --rm -it \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v ${PIPELINE_DIR}:/pipeline:ro \
    -v ${WORKDIR}:${WORKDIR} \
    -w ${WORKDIR} \
    docker.biotools.space/schemabio/nextflow:25.10.4 \
    nextflow run /pipeline/main.nf \
        -entry QC_ALIGNMENT \
        -config /pipeline/conf/qc_alignment.config \
        --config /pipeline/examples/qc_alignment.json \
        -profile docker
```

#### CNV 基线构建流程 CNV_BASELINE

使用正常样本构建 CNV 参考基线：

```bash
docker run --rm -it \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v ${PIPELINE_DIR}:/pipeline:ro \
    -v ${WORKDIR}:${WORKDIR} \
    -w ${WORKDIR} \
    docker.biotools.space/schemabio/nextflow:25.10.4 \
    nextflow run /pipeline/main.nf \
        -entry CNV_BASELINE \
        -config /pipeline/conf/cnv_baseline.config \
        --config /pipeline/examples/cnv_baseline.json \
        -profile docker
```

> **Docker-in-Docker 注意事项：** 除 pipeline 目录外，所有路径挂载时容器内外必须一致。因为 Nextflow 调度的子容器由宿主机 Docker daemon 启动，只能看到宿主机文件系统。

### 5. 配置文件说明

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

### 6. 输出结构

```
results/
└── {sample_id}/
    ├── 01_qc/              # 质控报告 (fastp, sex_check, coverage, metrics)
    ├── 02_alignment/       # 比对结果 (CRAM/CRAI)
    └── 03_variants/        # 变异结果
        ├── snv_indel/      # SNP/INDEL (DeepVariant, VEP, GenMod, Slivar)
        ├── mt/             # 线粒体变异 (Mutect2)
        ├── str/            # STR (ExpansionHunter, Stranger)
        └── cnv/            # CNV (CNVkit)
└── pipeline_info/          # 流程报告 (timeline/report/trace)
```

## 流程说明

### WES_SINGLE (默认完整流程)

```
FASTP (质控过滤) -> BWA_MEM2 (比对) -> GATK_MARKDUPLICATES (标记重复) -> CRAM
    -> DEEPVARIANT (SNV/INDEL) -> WHATSHAP (定相) -> VEP (注释) -> GENMOD (遗传模式)
    -> GATK_MUTECT2_MT (线粒体)
    -> EXPANSIONHUNTER (STR) -> STRANGER (STR注释)
    -> CNVKIT (CNV)
    -> PLINK2 (基因型分析) -> BCFTOOLS_ROH (ROH检测)
    -> SLIVAR (变异过滤)
```

### QC_ALIGNMENT (质控+比对)

```
FASTP (质控过滤) -> BWA_MEM2 (比对) -> CRAM
```

输出：质控报告 + sorted.cram

### CNV_BASELINE (CNV基线构建)

```
CNVKIT_TARGET + CNVKIT_ANTITARGET -> CNVKIT_REFERENCE_BUILD
```

输出：targets.bed + antitargets.bed + reference.cnn

## Profiles

| Profile | max_cpus | max_memory | 说明 |
|---------|----------|------------|------|
| `local` | 4 | 12 GB | 本地单机 (资源受限) |
| `docker` | - | - | Docker 容器运行时 |
| `singularity` | - | - | Singularity 容器运行时 |
| `podman` | - | - | Podman 容器运行时 |
| `slurm` | 32 | 128 GB | SLURM 调度器 |
| `lsf` | 32 | 128 GB | LSF 调度器 |
| `k8s` | 64 | 256 GB | Kubernetes |
| `aliyun` | 64 | 256 GB | 阿里云 |
| `tencent` | 64 | 256 GB | 腾讯云 |
| `test` | 2 | 6 GB | 测试用最小资源 |

可组合使用: `-profile local,docker` 或 `-profile slurm,singularity`

## 资源配置

资源通过命令行参数覆盖 profile 默认值：

```bash
# 限制最大资源
nextflow run /pipeline/main.nf \
    --config sample.json \
    --max_cpus 8 \
    --max_memory 16.GB \
    -profile docker
```

各模块的默认资源需求 (会被 max_* 参数限制)：

| 模块 | CPU | 内存 | 说明 |
|------|-----|------|------|
| FASTP | 8 | 4 GB | 多线程有效，内存需求低 |
| BWA_MEM2 | 16 | 64 GB | 内存 >= 64GB 用 bwa-mem2，否则自动降级为 bwa |
| GATK_MARKDUPLICATES | 4 | 16 GB | - |
| DEEPVARIANT | 16 | 32 GB | GPU 可加速 |
| VEP | 4 | 8 GB | - |

## 项目结构

```
schema-germline/
├── main.nf                 # 入口文件
├── nextflow.config         # 全局配置 (公共参数)
├── conf/
│   ├── modules.config      # 模块容器和发布配置
│   ├── wes_single.config   # WES_SINGLE 专用参数
│   ├── qc_alignment.config # QC_ALIGNMENT 专用参数
│   ├── cnv_baseline.config # CNV_BASELINE 专用参数
│   ├── aliyun.config       # 阿里云配置
│   └── tencent.config      # 腾讯云配置
├── modules/local/          # 本地模块
│   ├── fastp/main.nf       # 质控过滤
│   ├── bwa_mem2/main.nf    # 比对
│   └── gatk/main.nf        # GATK 工具集
├── workflows/
│   ├── wes_single.nf       # WES 单样本完整流程
│   ├── qc_alignment.nf    # 质控+比对流程
│   └── cnv_baseline.nf     # CNV 基线构建
├── examples/              # 示例配置文件
│   ├── sample.json         # WES_SINGLE 配置
│   ├── qc_alignment.json   # QC_ALIGNMENT 配置
│   └── cnv_baseline.json   # CNV_BASELINE 配置
└── containers/             # Dockerfile
```

## License

Apache License 2.0
