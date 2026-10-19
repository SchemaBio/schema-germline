# Schema Germline Pipeline

基于 WDL 的全外显子胚系变异检测流程，支持多种变异类型分析。

## 流程概览

| 流程 | 文件 | 用途 |
|------|------|------|
| SingleWES | `single.wdl` | 单样本全外显子分析 |
| VCF2Anno | `vcf2anno.wdl` | VCF 起始的 SNP/InDel 注释，支持 GRCh37/GRCh38 |
| CNVBaseline | `baseline.wdl` | CNV基线参考构建 |
| CNVBaselineFix | `baseline_fix.wdl` | 使用正常样本校正现有 CNV 基线 |
| Trio | `trio.wdl` | 家系父母子分析 |

`BwaAlign` 优先使用 minibwa：参考目录需包含与 FASTA 同名的 `.l2b`、`.mbw` 非空索引，运行节点需支持 AVX2/AVX512、具备 minibwa 且内存大于 16 GiB。选择不受线程数是否达到 16 限制，任务至少申请 32 GiB 内存。无法使用 minibwa 时，仅在标准 BWA 的 `.amb`、`.ann`、`.bwt`、`.pac`、`.sa` 索引完整且 `bwa` 可用时回退，否则明确报错。比对命令失败会使任务失败，不会被下游排序命令掩盖。

离线回归检查：`python -m unittest discover -s tests`（需要 Bash）。

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
docker pull docker.schema-bio.com/schemabio/germline:v0.1.3
docker pull docker.schema-bio.com/schemabio/gatk:4.6.2.0
docker pull docker.schema-bio.com/schemabio/deepvariant:1.10.0
docker pull docker.schema-bio.com/schemabio/vep:116.0
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

## VCF 起始的 SNP/InDel 注释

`vcf2anno.wdl` 提供 `VCF2Anno` 工作流，输入已有 VCF，复用当前 VEP 116 和 `GERMLINE.SNPInDelReport`。不需要 FASTQ、BAM、BED 或 CNV 基线，不重新检出变异或定相，不执行其他变异类型分析。输入 VCF 中的线粒体 SNP/InDel 可进入同一注释报表。

### 镜像与参数

沿用 `single.wdl` 的镜像内固定路径调用方式，新增两个脚本：

- `/opt/schema-germline/scripts/initialize_vcf.py`：样本提取、VCF 初始化和规范化。
- `/opt/schema-germline/scripts/infer_vcf_sex.py`：根据拆分前 X 染色体变异记录推算性别。

报表继续调用 `/opt/schema-germline/scripts/vep_report.py`。发布的 germline 镜像须封装当前版本的这三个脚本；不通过 WDL `File` 输入传入脚本。当前默认镜像名沿用 `docker.schema-bio.com/schemabio/germline:v0.1.7`，VEP 镜像沿用 `docker.schema-bio.com/schemabio/vep:116.0`。

hg19/hg38 使用同一份 WDL，不维护两套 JSON；参考配置由 Octopus 注入：

| 输入 | 说明 |
| --- | --- |
| `prefix` | 输出前缀，字母或数字开头，允许字母、数字、点、下划线、连字符，最多 128 字符 |
| `vcf` | 普通 VCF、gzip 或 BGZF 压缩 VCF，不要求输入索引 |
| `assembly` | hg19 对应 `GRCh37`；hg38 对应 `GRCh38` |
| `fasta` / `ref_dir` | 参考序列名称及其目录，必须存在对应 `.fai`；与 single 一样使用 `basename(fasta)` |
| `cache_dir` / `schema_bundle` | 对应组装的 VEP 116 merged cache 和现有注释库 |
| `sample_name` | 可选；单样本自动选择，多样本必须明确指定，不支持 sites-only VCF |
| `threads` | VEP worker 数，默认 8；VEP 内存沿用 `threads × 3 GiB` |
| `clinvar_version` | 默认 `20260415`，资源命名沿用当前 VEP 配置 |

参考库沿用当前流程的 ClinVar、cytoBand、gnomAD v4.1、Pangolin、EVOScore2 和 AlphaMissense v3 及所需索引。输入坐标、FASTA、cache 和资源必须属于同一组装，不执行 liftover；染色体名称统一不代表自动转换参考版本。

### 预处理与报表

初始化按 FASTA `.fai` 处理 `chr1/1`、`chrX/X`、`chrM/M/MT` 等别名，重建 contig 头信息。非标准 contig 必须与参考精确匹配，未知或歧义 contig、长度冲突、REF 不匹配均明确失败。

提取指定样本并去除未使用 ALT，随后使用 bcftools 进行 REF 校验、左对齐、共同碱基裁剪、多等位拆分、排序、BGZF 压缩和 TBI 索引。拆分重映射 GT/AD 等字段，保留源数据中可用的相位信息，不重新定相。原多等位描述记录在 `VCF2ANNO_ORIG` 中。

只保留该样本完整 GT 且实际携带的 SNP/InDel；跳过纯参考、缺失或部分缺失 GT、符号 ALT、`*`、参考区块及等长多碱基替换。一个记录内混合不同类型 ALT 时，只保留实际携带的受支持等位基因。多倍体或非法 GT 明确报错。不会将 MNV 原子化为 SNP，也不自动去重。

默认不按 FILTER 或 QUAL 丢弃记录，保留非 PASS、QUAL `.` 和 `0`。旧 CSQ 在重新注释前移除。VEP 使用自身并行，避免旧的染色体拆分/合并路径。注释发布前比较变异、GT、QUAL、FILTER，并核对报表数量，防止 VEP 或报表静默漏记录。

新流程对现有 `SNPInDelReport` 显式传入推算性别和 `external_vcf=true`，开启不过滤 QUAL、保留实际杂合 X/Y 的规则。原 single/trio 仍从 SRY 计数判断，默认过滤和合子状态行为不变。

### 性别推算口径

统计所选样本、拆分多等位前的有效 X 染色体 SNP/InDel 原始记录；每条输入记录最多计一次，包含 PAR 区域，不按 ALT 或转录本重复计数。

- `x_sum`：GT 完整且携带实际受支持 ALT 的 X 记录数。
- `x_hom`：二倍体纯合非参考记录数，如 `1/1`、`1|1`、`2/2`；`1/2` 不计入。
- `x_haploid_alt`：单倍体 ALT 记录数，如 `1`，单独记录。
- 比例为 `(x_hom + x_haploid_alt) / x_sum`。

默认 `sex_min_sites=100`、`female_max_ratio=0.5`、`male_min_ratio=0.7`。位点足够时，比例 ≤50% 推算女性，≥70% 推算男性，中间为未知；不足 100 个有效 X 位点或没有有效 X 位点也标记未知。三个参数均可配置，要求最小位点数 ≥1 且 `0 <= female_max_ratio < male_min_ratio <= 1`。

**由 VCF 的 X 染色体位点比例推算，可能受 X 染色体变异、检测区域及变异检出方式影响。** 这些经验阈值不等同于核型鉴定。结果保留计数、比例、阈值、推算结果、未知原因和该提示。

推算不修改 GT。实际杂合 X/Y 保留 `Het`；男性纯合 ALT 或显式单倍体 ALT 的 X/Y 可标记 `Hemi`。未知时二倍体纯合仍为 `Hom`，不默认套用女性或男性推断。

### 输出与验证

工作流输出规范化 VCF/TBI、VEP VCF/TBI、兼容现有列结构的 `.snv_indel.txt`、预处理统计 JSON、性别推算 JSON/TXT，以及包含各输出、统计与 `vep_skipped` 的 summary。没有有效变异时跳过 VEP，产生空 VCF、有效索引和仅有表头的报表，摘要明确标记跳过。

```bash
miniwdl check vcf2anno.wdl
python -m unittest discover -s tests
```

回归测试使用小型人工参考、真实 bcftools 和与 WDL 相同的命令验证初始化、推算及报表。非空报表测试使用合成 CSQ，不替代真实 hg19/hg38 注释库的完整 VEP 测试。测试应使用封装当前脚本的 germline 镜像；镜像外运行需安装 bcftools、pandas 和 openpyxl。

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
tar -zxvf homo_sapiens_merged_vep_116_GRCh38.tar.gz -C vep
rm homo_sapiens_merged_vep_116_GRCh38.tar.gz

# 可选，建立 bwa-mem2 索引
bwa-mem2 index Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

## License

MIT
