# Schema Germline Pipeline

基于 Nextflow DSL2 的全外显子胚系变异检测流程。

## 快速开始

### 1. 环境要求

- Docker
- Nextflow
- 自部署镜像已拉取（见下文）

### 2. 拉取镜像

```bash
# 流程依赖镜像
docker pull docker.schema-bio.com/schemabio/germline:v1.0.0
docker pull docker.schema-bio.com/schemabio/mapping:v1.0.0
docker pull docker.schema-bio.com/schemabio/gatk:4.6.2.0
docker pull docker.schema-bio.com/schemabio/cnvkit:0.9.13
docker pull docker.schema-bio.com/schemabio/expansionhunter:5.0.0
docker pull docker.schema-bio.com/schemabio/deepvariant:1.10.0
docker pull docker.schema-bio.com/schemabio/vep:115.2
docker pull docker.schema-bio.com/schemabio/glnexus:v1.4.1
docker pull docker.schema-bio.com/schemabio/whatshap:2.8
docker pull docker.schema-bio.com/schemabio/tiea_wes:2.0.0
docker pull docker.schema-bio.com/schemabio/automap:1.3
docker pull docker.schema-bio.com/schemabio/peddy:0.4.8
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

#### WES_SINGLE (单人完整流程)

```bash
# 设置路径
WORKDIR=/mnt/d/analysis
PIPELINE_DIR=/path/to/schema-germline

docker run --rm -it \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v ${PIPELINE_DIR}:/pipeline:ro \
    -v ${WORKDIR}:${WORKDIR} \
    -w ${WORKDIR} \
    docker.schema-bio.com/schemabio/nextflow:25.10.4 \
    nextflow run /pipeline/main.nf \
        --config /pipeline/examples/single.json \
        -profile local
```

#### WES_TRIO (家系 Trio 流程)

```bash
docker run --rm -it \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v ${PIPELINE_DIR}:/pipeline:ro \
    -v ${WORKDIR}:${WORKDIR} \
    -w ${WORKDIR} \
    docker.schema-bio.com/schemabio/nextflow:25.10.4 \
    nextflow run /pipeline/main.nf \
        -entry WES_TRIO \
        --config /pipeline/examples/trio.json \
        -profile local
```

#### CNV_BASELINE (CNV 基线构建)

```bash
docker run --rm -it \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v ${PIPELINE_DIR}:/pipeline:ro \
    -v ${WORKDIR}:${WORKDIR} \
    -w ${WORKDIR} \
    docker.schema-bio.com/schemabio/nextflow:25.10.4 \
    nextflow run /pipeline/main.nf \
        -entry CNV_BASELINE \
        --config /pipeline/examples/baseline.json \
        -profile local
```

> **说明：**
> - `-entry` 参数显式指定入口，不指定时自动检测配置文件类型
> - `--config` 参数指定 JSON 配置文件路径
> - `-profile` 参数指定运行环境

### 5. 配置文件说明

流程通过 JSON 配置文件指定输入参数，三种流程各有对应的配置格式：

#### WES_SINGLE 配置

```json
{
  "sample_id": "样本ID",
  "read1": "R1 fastq 路径",
  "read2": "R2 fastq 路径",
  "reference": {
    "fasta": "参考基因组 FASTA"
  },
  "outdir": "输出目录"
}
```

#### WES_TRIO 配置

```json
{
  "family_id": "家系ID",
  "samples": [
    {"sample_id": "先证者ID", "role": "proband", "read1": "...", "read2": "..."},
    {"sample_id": "父ID", "role": "father", "read1": "...", "read2": "..."},
    {"sample_id": "母ID", "role": "mother", "read1": "...", "read2": "..."}
  ],
  "reference": {"fasta": "..."},
  "outdir": "输出目录"
}
```

#### CNV_BASELINE 配置

```json
{
  "samples": [
    {"sample_id": "正常样本1", "read1": "...", "read2": "..."},
    {"sample_id": "正常样本2", "read1": "...", "read2": "..."}
  ],
  "reference": {"fasta": "..."},
  "target_bed": "捕获区域BED",
  "outdir": "输出目录"
}
```

> bwa/bwa-mem2 索引文件需与 fasta 同目录同前缀，无需单独配置。

## 流程说明

### WES_SINGLE (单人完整流程)

```
FASTP (质控过滤) -> BWA_MEM2 (比对) -> GATK_MARKDUPLICATES (标记重复) -> BAM
    -> DEEPVARIANT (SNV/INDEL) -> WHATSHAP (定相) -> VEP (注释)
    -> GATK_MUTECT2_MT (线粒体)
    -> EXPANSIONHUNTER (STR) -> STRANGER (STR注释)
    -> CNVKIT (CNV)
    -> TIEA_WES (MEI)
    -> AUTOMAP_ROH (ROH检测) -> AUTOMAP_UPD (UPD推断)
```

### WES_TRIO (家系 Trio 流程)

```
所有成员并行处理:
    FASTP (质控) -> BWA_MEM2 (比对) -> DEEPVARIANT (gVCF)

家系合并:
    GLNEXUS_JOINTCALL_TRIO (合并gVCF + de novo检测) -> PEDDY (亲缘验证)

合并VCF处理:
    WHATSHAP (家系定相) -> VEP (注释)

先证者专属分析:
    CNVKIT (CNV) + EXPANSIONHUNTER/STRANGER (STR) + TIEA_WES (MEI)
    + MUTECT2_MT (线粒体) + BAF (BAF矩阵)

Trio ROH/UPD:
    AUTOMAP_ROH (三人ROH) -> AUTOMAP_UPD (Trio精确UPD检测)
```

### CNV_BASELINE (CNV 基线构建)

```
BWA/BWA-MEM2 (比对) -> SAMTOOLS_INDEX (索引)
    -> CNVKIT_TARGET_ANTITARGET (生成BED)
    -> CNVKIT_COVERAGE (计算覆盖度)
    -> CNVKIT_REFERENCE_BUILD (构建基线)
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
| samtools | BAM 处理 | mapping:2026Jan |
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
| BCFTools | VCF 统计 | bcftools:1.21 |
| AutoMap | ROH 区域检测、近亲系数估算、UPD 推断 | automap:1.3 |

#### 过滤与筛选

| 软件 | 用途 | Docker 镜像 |
|------|------|-------------|

### Docker 镜像列表

所有镜像托管于 `docker.schema-bio.com/schemabio/`：

| 镜像 | 主要包含工具 |
|------|-------------|
| germline:v1.0.0 | 基础环境 |
| mapping:v1.0.0 | samtools, bwa, bwa-mem2, fastp, mosdepth, bamdst |
| gatk:4.6.2.0 | GATK + samtools |
| cnvkit:0.9.13 | CNVkit, Stranger, bedtools |
| expansionhunter:5.0.0 | ExpansionHunter |
| deepvariant:1.10.0 | DeepVariant |
| vep:115.2 | VEP |
| glnexus:v1.4.1 | GLnexus |
| whatshap:2.8 | WhatsHap |
| tiea_wes:2.0.0 | TIEA-WES (MEI检测) |
| automap:1.3 | AutoMap (ROH/UPD 分析) |
| peddy:0.4.8 | Peddy (亲缘验证) |

## Profiles

| Profile | max_cpus | max_memory | 说明 |
|---------|----------|------------|------|
| `local` | 16 | 64 GB | 本地单机 |
| `slurm` | 32 | 128 GB | SLURM 调度器 |
| `lsf` | 32 | 128 GB | LSF 调度器 |
| `tencent` | 64 | 256 GB | 腾讯云实例 |

Docker 容器配置见 `conf/docker.config`，定义了各模块使用的镜像。

### 运行环境配置

各运行环境的详细配置在独立的配置文件中：

| Profile | 配置文件 | 说明 |
|---------|---------|------|
| `local` | nextflow.config (内置) | 本地单机 |
| `slurm` | `conf/slurm.config` | SLURM HPC 集群 |
| `lsf` | `conf/lsf.config` | LSF HPC 集群 |
| `tencent` | `conf/tencent.config` | 腾讯云实例 |

修改示例 (`conf/slurm.config`)：

```groovy
process {
    executor       = 'slurm'
    queue          = 'your_queue'           // 修改为你的队列名称
    clusterOptions = '--account=your_account'  // 修改为你的账户
}

params {
    max_cpus   = 64                          // 根据集群资源调整
    max_memory = '256.GB'
    max_time   = '168.h'
}
```

### 自定义运行环境

1. **修改现有 Profile** - 直接编辑 `conf/slurm.config` 或 `conf/lsf.config`

2. **添加新 Profile** - 创建 `conf/my_cluster.config`：

```groovy
process {
    executor       = 'slurm'
    queue          = 'high_priority'
    clusterOptions = '--account=my_project'
}

params {
    max_cpus   = 64
    max_memory = '256.GB'
    max_time   = '168.h'
}
```

然后在 `nextflow.config` 的 `profiles` 块中添加：

```groovy
my_cluster {
    includeConfig 'conf/my_cluster.config'
}
```

运行时使用：`-profile my_cluster`

### 资源配置

资源可通过命令行参数覆盖 profile 默认值：

```bash
nextflow run main.nf \
    --config sample.json \
    --max_cpus 8 \
    --max_memory 16.GB \
    -profile local
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
├── main.nf                      # 主入口文件 (整合三个 workflow)
├── nextflow.config              # 全局配置
├── conf/                        # 配置文件
│   ├── modules.config           # 进程资源规格
│   ├── docker.config            # Docker 镜像配置
│   ├── slurm.config             # SLURM 集群配置
│   ├── lsf.config               # LSF 集群配置
│   ├── tencent.config           # 腾讯云实例配置
│   └── trio.config              # Trio 流程参数
├── modules/                    # 模块定义 (18个)
│   ├── fastp/main.nf           # 质控过滤
│   ├── bwamem/main.nf          # BWA/BWA-MEM2 比对
│   ├── samtools/main.nf        # BAM 处理
│   ├── sambamba/main.nf        # BAM 处理
│   ├── gatk/main.nf            # GATK 工具集
│   ├── deepvariant/main.nf     # SNV/INDEL 检测
│   ├── parabricks/main.nf      # GPU 加速比对和变异检测
│   ├── vep/main.nf             # 变异注释
│   ├── whatshap/main.nf        # 单倍型分相
│   ├── expansionhunter/main.nf # STR 检测
│   ├── stranger/main.nf        # STR 注释
│   ├── cnvkit/main.nf          # CNV 检测
│   ├── bcftools/main.nf        # VCF 工具
│   ├── tiea-wes/main.nf        # MEI 检测
│   ├── glnexus/main.nf         # gVCF 联合分型
│   ├── automap/main.nf         # ROH/UPD 分析
│   ├── peddy/main.nf           # 亲缘验证
│   └── xamdst/main.nf          # 覆盖度统计
├── workflows/                  # 工作流定义
│   ├── single/main.nf          # WES 单样本完整流程
│   ├── trio/main.nf            # WES Trio 家系流程
│   └── baseline/main.nf        # CNV 基线构建
├── examples/                    # 示例配置文件
│   ├── single.json              # WES_SINGLE 配置示例
│   ├── trio.json                # WES_TRIO 配置示例
│   └── baseline.json            # CNV_BASELINE 配置示例
└── README.md                    # 说明文档
```

### 输出结果目录结构

参考 Illumina DRAGEN、Sentieon、GATK Best Practices 的目录结构设计：
- 按数据类型分类，而非按步骤编号
- QC 报告集中管理，便于查看
- 中间结果 (raw) 和最终结果 (annotated) 分离
- 每种变异类型独立目录

```
results/
└── {sample_id}/
    ├── qc/                              # 质控报告 (集中管理)
    │   ├── fastp/                       # Fastp 质控
    │   │   ├── {sample}.json            # JSON 格式报告
    │   │   └── {sample}.html            # HTML 可视化报告
    │   ├── alignment/                   # 比对质控指标
    │   │   ├── {sample}.metrics.txt     # MarkDuplicates 指标
    │   │   └── {sample}_quality_metrics.*  # GATK CollectQCMetrics
    │   ├── coverage/                    # 覆盖度报告
    │   │   ├── {sample}.coverage.txt    # 常染色体覆盖度统计
    │   │   ├── {sample}_mt.coverage.txt # 线粒体覆盖度统计
    │   │   ├── {sample}.chromosome.txt  # 染色体覆盖度分布
    │   │   └── {sample}.region.txt      # 区域覆盖度统计
    │   └── sex_check/                   # 性别检测
    │       └── {sample}.sex.json        # 性别检测结果
    │
    ├── alignment/                       # 比对结果
    │   ├── {sample}.bam                 # 标记重复后的 BAM
    │   └── {sample}.bam.bai             # BAM 索引
    │
    ├── variants/                        # 变异检测结果
    │   ├── snv_indel/                   # SNV/Indel 变异
    │   │   ├── raw/                     # 原始检测结果
    │   │   │   ├── {sample}.vcf.gz      # DeepVariant VCF
    │   │   │   ├── {sample}.vcf.gz.tbi  # VCF 索引
    │   │   │   └── {sample}.g.vcf.gz    # gVCF (可选)
    │   │   ├── phased/                  # 单倍型定相结果
    │   │   │   ├── {sample}.phased.vcf.gz
    │   │   │   └── {sample}.phasing.json  # 定相统计
    │   │   └── annotated/               # VEP 注释结果 (最终)
    │   │       ├── {sample}.vep.vcf.gz
    │   │       ├── {sample}.vep.vcf.gz.tbi
    │   │       └── {sample}.vep_summary.txt
    │   │
    │   ├── cnv/                         # CNV 检测结果
    │   │   ├── {sample}.cnr             # CNVkit 拷贝数比率
    │   │   ├── {sample}.cns             # CNVkit 分段结果
    │   │   ├── {sample}.call.cns        # CNVkit 调用结果
    │   │   └── {sample}.cnv.vcf.gz      # CNV VCF 文件
    │   │
    │   ├── str/                         # STR 扩展检测
    │   │   ├── raw/                     # ExpansionHunter 原始输出
    │   │   │   ├── {sample}.vcf.gz
    │   │   │   └── {sample}.json        # STR profile
    │   │   ├── annotated/               # Stranger 注释结果
    │   │   │   ├── {sample}.stranger.vcf.gz
    │   │   │   └── {sample}.stranger_summary.json
    │   │   └ filtered/                 # 致病性 STR 过滤结果
    │   │       └── {sample}.pathogenic_str.vcf.gz
    │   │
    │   ├── mei/                         # MEI 检测结果
    │   │   ├── raw/                     # TIEA-WES 输出
    │   │   │   └── {sample}.mei.vcf.gz
    │   │   └ annotated/                 # VEP MEI 注释
    │   │       └── {sample}.mei.vep.vcf.gz
    │   │
    │   └── mt/                          # 线粒体变异
    │   │   ├── raw/                     # Mutect2 输出
    │   │   │   ├── {sample}.mt.vcf.gz
    │   │   │   └ {sample}.mt.stats
    │   │   └ annotated/                 # VEP MT 注释
    │   │       └── {sample}.mt.vep.vcf.gz
    │
    ├── baf/                             # BAF 矩阵 (辅助分析)
    │   ├── {sample}.baf.tsv             # BAF 矩阵 TSV
    │   └── {sample}.baf.json            # BAF 矩阵 JSON
    │
    ├── homozygosity/                    # ROH/UPD 分析结果
    │   ├── {sample}.roh.bed             # ROH 区域 BED
    │   ├── {sample}.roh.csv             # ROH 区域列表
    │   ├── {sample}.roh.json            # ROH 检测摘要
    │   ├── {sample}.roh.chromosome.txt  # 染色体分布报告
    │   ├── {sample}.inbreeding.json     # 近亲系数估算
    │   ├── {sample}.upd.bed             # UPD 区域 (推断)
    │   ├── {sample}.upd.json            # UPD 结果 JSON
    │   └── {sample}.upd_report.txt      # UPD 检测报告
    │
    └── pipeline_info/                   # Nextflow 流程信息
        ├── timeline.html                # 执行时间线
        ├── report.html                  # 执行报告
        ├── trace.txt                    # 追踪日志
        └── dag.svg                      # DAG 图
```

### CNV 基线目录结构 (CNV_BASELINE 流程)

```
results/
└── cnv_baseline/
    ├── targets.bed                      # 目标区域 BED
    ├── antitargets.bed                  # 反目标区域 BED
    └── reference.cnn                    # CNV 参考基线
```

### 家系 Trio 目录结构 (WES_TRIO 流程)

```
results/
└── {family_id}/
    ├── qc/                              # 所有成员质控报告
    │   ├── fastp/                       # Fastp 质控 (所有成员)
    │   ├── alignment/                   # 比对质控指标
    │   ├── coverage/                    # 覆盖度报告
    │   ├── sex_check/                   # 性别检测
    │   └── peddy/                       # 亲缘关系验证结果
    │
    ├── alignment/                       # 所有成员比对结果
    │   ├── {proband}.bam
    │   ├── {father}.bam
    │   ├── {mother}.bam
    │   └── *.bam.bai
    │
    ├── variants/                        # 家系变异检测结果
    │   ├── snv_indel/
    │   │   ├── raw/                     # 各成员 gVCF
    │   │   ├── joint/                   # GLnexus 合并 VCF + de novo VCF
    │   │   ├── phased/                  # WhatsHap 家系定相结果
    │   │   └── annotated/               # VEP 注释结果
    │   └── de_novo/                     # de novo 变异检测结果
    │
    ├── proband/                         # 先证者专属分析
    │   ├── cnv/                         # CNV 检测结果
    │   ├── str/                         # STR 检测结果
    │   │   ├── raw/
    │   │   ├── annotated/
    │   │   └ filtered/
    │   ├── mei/                         # MEI 检测结果
    │   │   ├── raw/
    │   │   └ annotated/
    │   ├── mt/                          # 线粒体变异
    │   │   ├── raw/
    │   │   └ annotated/
    │   └── baf/                         # BAF 矩阵
    │
    ├── homozygosity/                    # Trio ROH/UPD 分析
    │   ├── {proband}.roh.bed            # 先证者 ROH
    │   ├── {proband}.upd.bed            # 先证者 UPD (精确 Trio 检测)
    │   ├── {proband}.upd_report.txt     # UPD 报告
    │   └── *.inbreeding.json            # 近亲系数
    │
    └── pipeline_info/                   # Nextflow 流程信息
```

## License

Apache License 2.0