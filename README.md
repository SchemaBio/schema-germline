# Schema Germline Pipeline

基于 Nextflow DSL2 的胚系变异检测流程。

## 快速开始

### 1. 环境要求

- Docker / Podman / Singularity (根据运行环境选择)
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
docker pull docker.biotools.space/schemabio/core:v2.0.0
docker pull docker.biotools.space/schemabio/statistical:v2.0.0
docker pull docker.biotools.space/schemabio/base:v2.0.0
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

## 调用的软件

### 流程框架

| 软件 | 版本 | 用途 |
|------|------|------|
| Nextflow | DSL2 | 工作流引擎 |

### 容器运行时

- Docker / Podman / Singularity (根据 profile 选择)

### 生物信息学工具

#### 质控与比对

| 软件 | 用途 | Docker 镜像 |
|------|------|-------------|
| fastp | FASTQ 质控过滤 | mapping:2026Jan |
| bwa | 序列比对 (内存<64GB) | mapping:2026Jan |
| bwa-mem2 | 序列比对 (内存≥64GB) | mapping:2026Jan |
| samtools | BAM/CRAM 处理 | mapping:2026Jan |
| bamdst | 覆盖度统计 | mapping:2026Jan |

#### 变异检测

| 软件 | 用途 | Docker 镜像 |
|------|------|-------------|
| DeepVariant | SNV/INDEL 检测 | deepvariant:1.10.0 |
| GATK Mutect2 | 线粒体变异检测 | gatk:4.6.2.0 |
| GATK MarkDuplicates | 重复标记 | gatk:4.6.2.0 |
| GATK CollectQCMetrics | 质控指标 | gatk:4.6.2.0 |
| ExpansionHunter | STR 检测 | expansionhunter:5.0.0 |
| CNVkit | CNV 检测与基线构建 | cnvkit:v0.9.11.p4 |
| GLnexus | gVCF 联合分型 | glnexus:v1.4.1 |

#### 定相与注释

| 软件 | 用途 | Docker 镜像 |
|------|------|-------------|
| WhatsHap | 单倍型分相 | core:v2.0.0 |
| VEP (Ensembl) | 变异注释 | ensembl-vep:release_115.2 |
| Stranger | STR 注释 | cnvkit:v0.9.11.p4 |
| SVDB | CNV/SV 数据库注释 | cnvkit:v0.9.11.p4 |
| GenMod | 遗传模式分析 | core:v2.0.0 |

#### 统计分析

| 软件 | 用途 | Docker 镜像 |
|------|------|-------------|
| Plink2 | 基因型分析、PCA、亲缘关系 | statistical:v2.0.0 |
| BCFTools | ROH 检测、VCF 统计 | statistical:v2.0.0 |

#### 过滤与筛选

| 软件 | 用途 | Docker 镜像 |
|------|------|-------------|
| Slivar | 变异过滤、复合杂合检测 | core:v2.0.0 |

#### 辅助脚本

| 软件 | 用途 | Docker 镜像 |
|------|------|-------------|
| aggregate_qc.py | QC 指标汇总 | base:v2.0.0 |
| results_to_parquet.py | 结果转 Parquet | base:v2.0.0 |
| check_fingerprint.py | 样本指纹验证 | base:v2.0.0 |
| make_pedigree.py | 家系文件生成 | base:v2.0.0 |

### Docker 镜像列表

所有镜像托管于 `docker.biotools.space/schemabio/`：

| 镜像 | 主要包含工具 |
|------|-------------|
| nextflow:25.10.4 | Nextflow 运行环境 |
| mapping:2026Jan | samtools, bwa, bwa-mem2, fastp, mosdepth, bamdst |
| gatk:4.6.2.0 | GATK + samtools |
| cnvkit:v0.9.11.p4 | CNVkit, Stranger, SVDB, bedtools |
| expansionhunter:5.0.0 | ExpansionHunter, Surf, Manta |
| deepvariant:1.10.0 | DeepVariant |
| ensembl-vep:release_115.2 | VEP |
| glnexus:v1.4.1 | GLnexus |
| core:v2.0.0 | GenMod, Slivar, WhatsHap |
| statistical:v2.0.0 | Plink2, Vcftools |
| base:v2.0.0 | Python 基础环境 (辅助脚本) |

## Profiles

| Profile | max_cpus | max_memory | 说明 |
|---------|----------|------------|------|
| `local` | 4 | 12 GB | 本地单机 (资源受限) |
| `docker` | 16 | 64 GB | Docker 容器运行时 |
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

## 目录结构

### 项目源码结构

```
schema-germline/
├── main.nf                      # 入口文件 (3个工作流入口)
├── nextflow.config              # 全局配置
├── conf/                        # 配置文件
│   ├── modules.config           # 模块容器和发布配置
│   ├── wes_single.config        # WES_SINGLE 专用参数
│   ├── qc_alignment.config      # QC_ALIGNMENT 专用参数
│   ├── cnv_baseline.config      # CNV_BASELINE 专用参数
│   ├── aliyun.config            # 阿里云配置
│   └── tencent.config           # 腾讯云配置
├── modules/local/               # 本地模块 (19个)
│   ├── fastp/main.nf            # 质控过滤
│   ├── bwa_mem2/main.nf         # 比对
│   ├── samtools/main.nf         # BAM/CRAM 处理
│   ├── gatk/main.nf             # GATK 工具集
│   ├── deepvariant/main.nf      # SNV/INDEL 检测
│   ├── vep/main.nf              # 变异注释
│   ├── genmod/main.nf           # 遗传模式分析
│   ├── whatshap/main.nf         # 单倍型分相
│   ├── expansionhunter/main.nf  # STR 检测
│   ├── stranger/main.nf         # STR 注释
│   ├── cnvkit/main.nf           # CNV 检测
│   ├── svdb/main.nf             # CNV/SV 注释
│   ├── plink2/main.nf           # 基因型分析
│   ├── bcftools/main.nf         # ROH/VCF 统计
│   ├── slivar/main.nf           # 变异过滤
│   ├── glnexus/main.nf          # gVCF 联合分型
│   ├── sample_identity/main.nf  # 样本指纹验证
│   ├── sex_check/main.nf        # 性别检查
│   └── xamdst/main.nf           # 覆盖度统计
├── workflows/                   # 工作流定义
│   ├── wes_single.nf            # WES 单样本完整流程
│   ├── qc_alignment.nf          # 质控+比对流程
│   └── cnv_baseline.nf          # CNV 基线构建
├── containers/                  # Dockerfile (18个)
│   ├── base/Dockerfile          # Python 基础环境
│   ├── mapping/Dockerfile       # 比对工具集
│   ├── gatk/Dockerfile          # GATK 工具集
│   ├── cnvkit/Dockerfile        # CNV 分析工具集
│   ├── deepvariant/Dockerfile   # DeepVariant
│   ├── vep/Dockerfile           # VEP
│   ├── expansionhunter/Dockerfile # STR 检测工具集
│   ├── glnexus/Dockerfile       # GLnexus
│   ├── core/Dockerfile          # 核心分析工具集
│   ├── statistical/Dockerfile   # 统计工具集
│   └── ... (其他)
├── bin/                         # Python 辅助脚本
│   ├── aggregate_qc.py          # QC 指标汇总
│   ├── bam_cram_converter.py    # BAM/CRAM 转换
│   ├── cnv_segment.py           # CNV 分段处理
│   ├── make_pedigree.py         # 家系文件生成
│   ├── results_to_parquet.py    # 结果转 Parquet
│   └── check_fingerprint.py     # 样本指纹验证
├── examples/                    # 示例配置文件
│   ├── sample.json              # WES_SINGLE 配置示例
│   ├── qc_alignment.json        # QC_ALIGNMENT 配置示例
│   └── cnv_baseline.json        # CNV_BASELINE 配置示例
├── docs/                        # 文档
├── assets/                      # 资源文件
├── design.md                    # 设计文档
└── README.md                    # 说明文档
```

### 输出结果目录结构

```
results/
└── {sample_id}/
    ├── 01_qc/                   # 质控报告
    │   ├── fastp.json           # fastp 质控 JSON
    │   ├── fastp.html           # fastp 质控 HTML 报告
    │   ├── *.metrics.txt        # GATK MarkDuplicates 指标
    │   ├── *.coverage*.txt      # bamdst 覆盖度统计
    │   ├── *.chromosome*.txt    # 染色体覆盖度
    │   ├── *.region*.txt        # 区域覆盖度
    │   └── *.sex*.txt           # 性别检查结果
    ├── 02_alignment/            # 比对结果
    │   ├── *.marked.cram        # 去重后的 CRAM
    │   └── *.marked.cram.crai   # CRAM 索引
    ├── 03_variants/             # 变异结果
    │   ├── snv_indel/           # SNP/INDEL
    │   │   ├── *.vcf.gz         # DeepVariant VCF
    │   │   ├── *.g.vcf.gz       # gVCF (可选)
    │   │   ├── *.phase.vcf.gz   # WhatsHap 定相结果
    │   │   ├── *.vep.vcf.gz     # VEP 注释结果
    │   │   ├── *.genmod*.vcf.gz # GenMod 遗传模式
    │   │   ├── *.slivar.vcf.gz  # Slivar 过滤结果
    │   │   ├── *.compound.vcf.gz # 复合杂合
    │   │   ├── *.joint.vcf.gz   # GLnexus 联合分型
    │   │   ├── *.pgen           # Plink2 pgen
    │   │   ├── *.pvar           # Plink2 pvar
    │   │   ├── *.psam           # Plink2 psam
    │   │   ├── *.king*          # 亲缘关系
    │   │   ├── *.eigenvec       # PCA 结果
    │   │   ├── *.roh*           # ROH 检测
    │   │   └── *.vcfstats.txt   # VCF 统计
    │   ├── mt/                  # 线粒体变异
    │   │   ├── *.vcf.gz         # Mutect2 VCF
    │   │   └── *.genmod*.vcf.gz # GenMod 注释
    │   ├── str/                 # STR
    │   │   ├── *.vcf.gz         # ExpansionHunter VCF
    │   │   ├── *.stranger.vcf.gz # Stranger 注释
    │   │   └── *.json           # STR profile
    │   └── cnv/                 # CNV
    │       ├── *.cnr            # CNVkit 拷贝数比率
    │       ├── *.cns            # CNVkit 分段
    │       ├── *.call.cns       # CNVkit 调用分段
    │       ├── *.vcf.gz         # CNV VCF
    │       └── *.svdb.vcf.gz    # SVDB 注释
    └── 04_reports/              # 输出报告
        ├── *.qc_summary.*       # QC 汇总
        ├── *.slivar.tsv.gz      # 变异表格
        ├── *.parquet            # Parquet 格式结果
        └── *.duodel.tsv.gz      # 复合杂合表格
└── cnv_baseline/                # CNV 基线 (仅 CNV_BASELINE)
│   ├── targets.bed              # 目标区域
│   ├── antitargets.bed          # 反目标区域
│   └── reference.cnn            # CNV 参基线
└── pipeline_info/               # 流程执行报告
    ├── timeline.html            # 时间线
    ├── report.html              # 执行报告
    ├── trace.txt                # 追踪日志
    └── dag.svg                  # DAG 图
```

## License

Apache License 2.0