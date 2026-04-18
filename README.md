# Schema Germline Pipeline

基于 Nextflow DSL2 的全外显子胚系变异检测流程。

## 快速开始

### 1. 环境要求

- Docker
- Nextflow
- 自部署镜像已拉取（见下文）

### 2. 拉取镜像

```bash
# 流程依赖的其他镜像（自动按需拉取）
docker pull docker.schema-bio.com/schemabio/mapping:v1.0.0
docker pull docker.schema-bio.com/schemabio/gatk:4.6.2.0
docker pull docker.schema-bio.com/schemabio/cnvkit:v0.9.13
docker pull docker.schema-bio.com/schemabio/expansionhunter:5.0.0
docker pull docker.schema-bio.com/schemabio/deepvariant:1.10.0
docker pull docker.schema-bio.com/schemabio/vep:115.2
docker pull docker.schema-bio.com/schemabio/glnexus:v1.4.1
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
    docker.schema-bio.com/schemabio/nextflow:25.10.4 \
    nextflow run /pipeline/main.nf \
        -config /pipeline/conf/wes_single.config \
        --config /pipeline/examples/sample.json \
        -profile docker
```

#### 质控+比对流程 QC_ALIGNMENT

仅进行 FASTP 质控和 BWA 比对，输出 sorted.bam：

```bash
docker run --rm -it \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v ${PIPELINE_DIR}:/pipeline:ro \
    -v ${WORKDIR}:${WORKDIR} \
    -w ${WORKDIR} \
    docker.schema-bio.com/schemabio/nextflow:25.10.4 \
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
    docker.schema-bio.com/schemabio/nextflow:25.10.4 \
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

### WES_SINGLE (默认完整流程 - 单人)

```
FASTP (质控过滤) -> BWA_MEM2 (比对) -> GATK_MARKDUPLICATES (标记重复) -> BAM
    -> DEEPVARIANT (SNV/INDEL) -> WHATSHAP (定相) -> VEP (注释)
    -> GATK_MUTECT2_MT (线粒体)
    -> EXPANSIONHUNTER (STR) -> STRANGER (STR注释)
    -> CNVKIT (CNV)
    -> TIEA_WES (MEI)
    -> AUTOMAP_ROH (ROH检测) -> AUTOMAP_UPD (UPD推断)
```

### WES_TRIO (家系流程 - 父母子 Trio)

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

### QC_ALIGNMENT (质控+比对)

```
FASTP (质控过滤) -> BWA_MEM2 (比对) -> BAM
```

输出：质控报告 + sorted.bam

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
| Plink2 | 基因型分析、PCA、亲缘关系 | statistical:v2.0.0 |
| BCFTools | ROH 检测、VCF 统计 | statistical:v2.0.0 |
| AutoMap | ROH 区域检测、近亲系数估算、UPD 推断 | automap:1.3 |

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

所有镜像托管于 `docker.schema-bio.com/schemabio/`：

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
| automap:1.3 | AutoMap (ROH/UPD 分析) |
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
│   ├── samtools/main.nf         # BAM 处理
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