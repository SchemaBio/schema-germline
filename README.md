# Schema Germline Pipeline

基于 WDL 的全外显子胚系变异检测流程，支持多种变异类型分析。

## 流程概览

| 流程 | 文件 | 用途 |
|------|------|------|
| SingleWES | `single.wdl` | 单样本全外显子分析 |
| CNVBaseline | `baseline.wdl` | CNV基线参考构建 |
| Trio | `trio.wdl` | 家系三人组分析（开发中） |

### SingleWES 分析模块

| 模块 | 功能 |
|------|------|
| Fastp | FASTQ质控过滤 |
| BWA-MEM | 序列比对 |
| Sambamba | 标记重复 |
| DeepVariant | SNP/InDel检测 |
| WhatsHap | 单倍型定相 |
| VEP | 变异注释（含ClinVar） |
| GATK Mutect2 | 线粒体变异检测 |
| CNVKit | 拷贝数变异分析 |
| ExpansionHunter | STR扩增检测 |
| Stranger | STR致病性注释 |
| AutoMap | ROH/UPD分析 |
| TIEA-WES | 转座子插入分析 |
| GATK/QCMetrics | 覆盖度质控指标 |

## 快速开始

### 1. 环境要求

- Docker
- miniwdl (`pip install miniwdl`)
- 参考基因组数据库

#### HPC 集群
- Slurm (`pip install miniwdl-slurm`)
- LSF (`pip install miniwdl-lsf`)

### 2. 拉取镜像

```bash
docker pull docker.schema-bio.com/schemabio/gatk:4.6.2.0
docker pull docker.schema-bio.com/schemabio/deepvariant:1.10.0
docker pull docker.schema-bio.com/schemabio/vep:115.2
docker pull docker.schema-bio.com/schemabio/whatshap:2.8
docker pull docker.schema-bio.com/schemabio/mapping:v1.0.0
docker pull docker.schema-bio.com/schemabio/germline:v0.0.4
```

### 3. 准备配置文件

#### 单样本分析 (`inputs/single.json`)

```json
{
    "SingleWES.prefix": "proband",
    "SingleWES.read_1": "/mnt/data/test/sample_R1.fq.gz",
    "SingleWES.read_2": "/mnt/data/test/sample_R2.fq.gz",
    "SingleWES.fasta": "/database/Homo_sapiens.GRCh37.dna.primary_assembly.fa",
    "SingleWES.bed": "/mnt/data/test/capture.bed",
    "SingleWES.flank_size": 50,
    "SingleWES.assembly": "GRCh37",
    "SingleWES.ref_dir": "/database",
    "SingleWES.cnvkit_reference": "/database/reference.cnvkit.cnn",
    "SingleWES.cnvkit_segmentation_method": "hmm-germline",
    "SingleWES.sry_sex_cutoff": 20,
    "SingleWES.cache_dir": "/database/vep",
    "SingleWES.schema_bundle": "/database/schema_bundle"
}
```

#### CNV基线构建 (`inputs/baseline.json`)

```json
{
    "CNVBaseline.prefix": "reference",
    "CNVBaseline.bed": "/mnt/data/test/capture.bed",
    "CNVBaseline.fasta": "/database/Homo_sapiens.GRCh37.dna.primary_assembly.fa",
    "CNVBaseline.assembly": "GRCh37",
    "CNVBaseline.read_1": [
        "/mnt/data/test/control1_R1.fq.gz",
        "/mnt/data/test/control2_R1.fq.gz"
    ],
    "CNVBaseline.read_2": [
        "/mnt/data/test/control1_R2.fq.gz",
        "/mnt/data/test/control2_R2.fq.gz"
    ],
    "CNVBaseline.ref_dir": "/database"
}
```

> **注意**: BWA索引文件需与fasta同目录同前缀，无需单独配置。

### 4. 配置 miniwdl

| 配置文件 | 环境 |
|----------|------|
| `conf/local.cfg` | 本地单机 |
| `conf/slurm.cfg` | Slurm集群 |
| `conf/lsf.cfg` | LSF集群 |

主要配置项：
- `max_tasks`: 最大并行任务数
- `defaults`: 默认Docker镜像和资源
- `root`: 工作目录路径

### 5. 运行流程

#### 单样本分析

```bash
miniwdl run single.wdl \
    --cfg conf/local.cfg \
    -i inputs/single.json \
    --dir /mnt/data/output \
    -p /path/to/schema-germline
```

#### CNV基线构建

```bash
miniwdl run baseline.wdl \
    --cfg conf/local.cfg \
    -i inputs/baseline.json \
    --dir /mnt/data/output \
    -p /path/to/schema-germline
```

#### HPC集群运行

```bash
# Slurm
miniwdl run single.wdl --cfg conf/slurm.cfg -i inputs/single.json

# LSF
miniwdl run single.wdl --cfg conf/lsf.cfg -i inputs/single.json
```

## 项目结构

```
schema-germline/
├── single.wdl          # 单样本流程
├── baseline.wdl        # CNV基线流程
├── trio.wdl            # 家系流程
├── tasks/              # 任务模块
│   ├── fastp.wdl       # 质控过滤
│   ├── bwamem.wdl      # BWA比对
│   ├── sambamba.wdl    # 标记重复
│   ├── samtools.wdl    # SAMtools工具
│   ├── xamdst.wdl      # 质控统计
│   ├── gatk.wdl        # GATK工具
│   ├── deepvariant.wdl # DeepVariant
│   ├── whatshap.wdl    # 定相
│   ├── vep.wdl         # 变异注释
│   ├── cnvkit.wdl      # CNV分析
│   ├── expansionhunter.wdl # STR检测
│   ├── stranger.wdl    # STR注释
│   ├── automap.wdl     # ROH/UPD
│   ├── tiea_wes.wdl    # 转座子分析
│   ├── germline.wdl    # 通用工具
│   ├── peddy.wdl       # 家系验证
│   ├── glnexus.wdl     # 变异合并
│   └── result.wdl      # 结果整理
├── conf/               # 配置文件
│   ├── local.cfg
│   ├── slurm.cfg
│   └── lsf.cfg
├── inputs/             # 输入示例
│   ├── single.json
│   └── baseline.json
└── assets/             # 资源文件
    └── mito.bed
```

## 参考数据库要求

- 参考基因组 (GRCh37/GRCh38)
- BWA索引文件
- VEP缓存目录
- ClinVar数据库
- CNVKit参考文件（可选，用于CNV分析）
- STR目录文件（ExpansionHunter）
- Schema Bundle（VEP插件资源）

## License

MIT