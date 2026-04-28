# Schema Germline Pipeline

基于 WDL 的全外显子胚系变异检测流程，支持多种变异类型分析。

## 流程概览

| 流程 | 文件 | 用途 |
|------|------|------|
| SingleWES | `single.wdl` | 单样本全外显子分析 |
| CNVBaseline | `baseline.wdl` | CNV基线参考构建 |
| Trio | `trio.wdl` | 家系父母子分析 |

## 快速开始

### 1. 环境要求

- Docker
- miniwdl (`pip install miniwdl`)
- 参考基因组数据库

拉取流程
```bash
git clone https://github.com/schemabio/schema-germline.git
```

#### HPC 集群
- Slurm (`pip install miniwdl-slurm`)
- LSF (`pip install miniwdl-lsf`)


#### **参数调整**

建议自行调整`single.wdl`、`trio.wdl`、`baseline.wdl`、`conf/*.cfg`的配置细则，适配自己的环境。

### 2. 拉取镜像

```bash
# 核心镜像
docker pull docker.schema-bio.com/schemabio/germline:v0.1.1
docker pull docker.schema-bio.com/schemabio/gatk:4.6.2.0
docker pull docker.schema-bio.com/schemabio/deepvariant:1.10.0
docker pull docker.schema-bio.com/schemabio/deeptrio:1.10.0
docker pull docker.schema-bio.com/schemabio/vep:115.2
docker pull docker.schema-bio.com/schemabio/whatshap:2.8
docker pull docker.schema-bio.com/schemabio/mapping:v1.0.0
docker pull docker.schema-bio.com/schemabio/cnvkit:0.9.13.2
docker pull docker.schema-bio.com/schemabio/cnvanno:v0.0.2
docker pull docker.schema-bio.com/schemabio/expansionhunter:5.0.0
docker pull docker.schema-bio.com/schemabio/stranger:v0.10.0.1
docker pull docker.schema-bio.com/schemabio/automap:1.3
docker pull docker.schema-bio.com/schemabio/tiea_wes:2.0.1
```

### 3. 准备配置文件

输入参数配置文件示例位于 `inputs/` 目录：

| 文件 | 用途 |
|------|------|
| `inputs/single.json` | 单样本分析参数 |
| `inputs/baseline.json` | CNV基线构建参数 |
| `inputs/trio.json` | 家系父母子构建参数 |

miniwdl运行配置文件示例位于 `conf/` 目录：

| 文件 | 环境 |
|------|------|
| `conf/local.cfg` | 本地单机 |
| `conf/slurm.cfg` | Slurm集群 |
| `conf/lsf.cfg` | LSF集群 |

> **注意**: BWA索引文件需与fasta同目录同前缀，无需单独配置。

### 4. 运行流程

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

#### 家系父母子分析

```bash
miniwdl run trio.wdl \
    --cfg conf/local.cfg \
    -i inputs/trio.json \
    --dir /mnt/data/output \
    -p /path/to/schema-germline
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

## 参考数据库下载

按需在 [HuggingFace](https://huggingface.co/datasets/pzweuj/SchemaBio_Bundle) 下载对应版本的参考数据库，以 hg38 为例：

```bash
# 使用 huggingface-cli 下载（推荐）
pip install huggingface_hub
huggingface-cli download pzweuj/SchemaBio_Bundle --repo-type dataset --include "hg38/*" --local-dir ./database

# 或使用 git clone
git clone https://huggingface.co/datasets/pzweuj/SchemaBio_Bundle ./database
```

建议储存路径：

```bash
cd database

# 解压参考基因组
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Schema Bundle
mkdir -p schema_bundle
mv hg38* schema_bundle/

# Bed
mkdir -p bed
mv Gencode.GRCh38.cnvkit.target.bed bed/

# 解压 VEP cache
mkdir vep
tar -zxvf homo_sapiens_merged_vep_115_GRCh38.tar.gz -C vep
rm homo_sapiens_merged_vep_115_GRCh38.tar.gz

# 可选，建立 bwa-mem2 索引
bwa-mem2 index Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

## License

MIT