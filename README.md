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
docker pull docker.schema-bio.com/schemabio/tiea_wes:2.0.1
docker pull docker.schema-bio.com/schemabio/automap:1.3
docker pull docker.schema-bio.com/schemabio/peddy:0.4.8
```

### 3. 准备配置文件

每个样本一个 JSON 文件，如 `examples/single.json`:

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

```bash
nextflow run main.nf \
    --config examples/single.json \
    -profile local
```

## WES_SINGLE 单人流程分析线路图

```mermaid
flowchart TB
    subgraph Input["输入数据"]
        FASTQ["FASTQ (R1/R2)"]
        REF["参考基因组 (FASTA)"]
        BED["捕获区域 BED"]
        SNP_POS["SNP 位点文件"]
        CNV_REF["CNVkit 基线"]
        STR_CAT["STR Catalog"]
    end

    subgraph QC["Step 1: 质控过滤"]
        FASTP["FastP\n接头去除 + 质控"]
    end

    subgraph Alignment["Step 2-3: 比对与后处理"]
        direction TB
        AlignBranch{"资源判断"}
        
        subgraph GPU_Mode["GPU 模式"]
            PB_FQ2BAM["Parabricks FQ2BAM\n(自带 MarkDup)"]
        end
        
        subgraph High_Mode["高资源模式 (≥64GB)"]
            BWAMEM2["BWA-MEM2"]
            MARKDUP1["GATK MarkDuplicates"]
            IDX1["Samtools Index"]
        end
        
        subgraph Low_Mode["低资源模式"]
            BWAMEM["BWA-MEM"]
            MARKDUP2["GATK MarkDuplicates"]
            IDX2["Samtools Index"]
        end
        
        BAM["BAM 文件"]
    end

    subgraph QC_Stats["Step 4-7: QC 统计"]
        direction TB
        XAMDST["XAMDST\n覆盖度统计"]
        XAMDST_MT["XAMDST (MT)\n线粒体覆盖度"]
        COLLECTQC["GATK CollectMultipleMetrics"]
        SEX_CHECK["性别检测 (SRY)"]
        BAF["BAF 矩阵计算"]
    end

    subgraph VariantCalling["Step 8: 变异检测"]
        direction TB
        DVBranch{"GPU 判断"}
        
        subgraph GPU_DV["GPU 模式"]
            PB_DV["Parabricks DeepVariant"]
        end
        
        subgraph CPU_DV["CPU 模式"]
            DV["DeepVariant (WES)"]
        end
        
        VCF["VCF/gVCF"]
        MT_MT2["GATK Mutect2 (MT)\n线粒体变异"]
    end

    subgraph Phasing["Step 9: 单倍型定相"]
        WHATSHAP["WhatsHap Phase"]
        PHASED_VCF["定相 VCF"]
    end

    subgraph SpecialVariants["Step 10-12: 特殊变异检测"]
        direction TB
        
        subgraph CNV_Detection["CNV 检测"]
            CNVKIT["CNVkit Batch"]
            CNV_VCF["CNV VCF"]
        end
        
        subgraph STR_Detection["STR 扩展检测"]
            EH["ExpansionHunter"]
            STRANGER["Stranger 注释"]
            STR_FILTER["致病性筛选"]
            STR_VCF["STR VCF"]
        end
        
        subgraph MEI_Detection["MEI 检测"]
            TIEA["TIEA-WES"]
            MEI_VCF["MEI VCF"]
        end
    end

    subgraph Annotation["Step 14: 变异注释"]
        direction TB
        
        subgraph VEP_Annotation["VEP 注释"]
            VEP_SNV["VEP (SNV/Indel)\ngnomAD/AlphaMissense\nClinVar/EVOScore"]
            VEP_MT["VEP (MT)\n线粒体注释"]
            VEP_MEI["VEP (MEI)\n移动元件注释"]
        end
    end

    subgraph Homozygosity["Step 15-16: 纯合性分析"]
        direction TB
        ROH["AutoMap ROH\n纯合区域检测"]
        UPD["AutoMap UPD\n单亲二倍体推断"]
    end

    subgraph Output["输出结果"]
        direction TB
        OUT_QC["QC 报告\n(FastP/GATK/覆盖度)"]
        OUT_VCF["变异 VCF\n(SNV/CNV/STR/MEI/MT)"]
        OUT_ANNOTATED["注释 VCF\n(VEP annotated)"]
        OUT_ROH["ROH/UPD 报告\n(近亲系数)"]
    end

    %% 输入连接
    FASTQ --> FASTP
    REF --> PB_FQ2BAM & BWAMEM2 & BWAMEM & DV & PB_DV
    BED --> XAMDST & DV & PB_DV & CNVKIT & ROH
    SNP_POS --> BAF
    CNV_REF --> CNVKIT
    STR_CAT --> EH & STRANGER

    %% 质控流程
    FASTP --> AlignBranch
    
    %% 比对分支
    AlignBranch -->|"use_gpu=true"| GPU_Mode
    AlignBranch -->|"use_gpu=false\n≥64GB"| High_Mode
    AlignBranch -->|"use_gpu=false\n<64GB"| Low_Mode
    
    PB_FQ2BAM --> BAM
    BWAMEM2 --> MARKDUP1 --> IDX1 --> BAM
    BWAMEM --> MARKDUP2 --> IDX2 --> BAM

    %% QC 统计 (并行)
    BAM --> XAMDST & XAMDST_MT & COLLECTQC & SEX_CHECK & BAF

    %% 变异检测
    BAM --> DVBranch
    DVBranch -->|"use_gpu=true"| GPU_DV
    DVBranch -->|"use_gpu=false"| CPU_DV
    PB_DV --> VCF
    DV --> VCF
    BAM --> MT_MT2

    %% 定相
    VCF --> WHATSHAP --> PHASED_VCF

    %% 特殊变异检测 (并行)
    BAM --> CNVKIT --> CNV_VCF
    BAM --> EH --> STRANGER --> STR_FILTER --> STR_VCF
    BAM --> TIEA --> MEI_VCF

    %% VEP 注释
    PHASED_VCF --> VEP_SNV
    MT_MT2 --> VEP_MT
    MEI_VCF --> VEP_MEI

    %% 纯合性分析
    VEP_SNV --> ROH & UPD

    %% 输出
    FASTP --> OUT_QC
    XAMDST & COLLECTQC --> OUT_QC
    VCF & CNV_VCF & STR_VCF & MEI_VCF --> OUT_VCF
    VEP_SNV & VEP_MT & VEP_MEI --> OUT_ANNOTATED
    ROH & UPD --> OUT_ROH
```

### 流程步骤详解

| 步骤 | 模块 | 功能 | 输出 |
|------|------|------|------|
| 1 | **FastP** | 接头去除、质控过滤 | clean FASTQ, QC报告 |
| 2-3 | **BWA/BWA-MEM2/Parabricks** | 序列比对 + MarkDuplicates | BAM文件 |
| 4 | **XAMDST** | 覆盖度统计 | coverage report |
| 5 | **GATK CollectMultipleMetrics** | 比对质量统计 | QC metrics |
| 6 | **SEX_CHECK_SRY** | 性别检测 | sex JSON |
| 7 | **BCFTOOLS BAF** | BAF矩阵计算 | BAF TSV/JSON |
| 8 | **DeepVariant** | SNV/Indel检测 | VCF/gVCF |
| 8b | **GATK Mutect2 (MT)** | 线粒体变异检测 | MT VCF |
| 9 | **WhatsHap** | 单倍型定相 | phased VCF |
| 10 | **CNVkit** | CNV检测 | CNV VCF |
| 11 | **ExpansionHunter + Stranger** | STR扩展检测 + 注释 | STR VCF |
| 12 | **TIEA-WES** | MEI检测 | MEI VCF |
| 14 | **VEP** | 变异功能注释 | annotated VCF |
| 15 | **AutoMap ROH** | 纯合区域检测 | ROH BED, 近亲系数 |
| 16 | **AutoMap UPD** | 单亲二倍体推断 | UPD report |

### 资源模式选择

流程会根据配置自动选择最优比对策略：

| 条件 | 比对模式 | 变异检测模式 |
|------|---------|-------------|
| `use_gpu=true` | Parabricks FQ2BAM | Parabricks DeepVariant |
| `use_gpu=false` + 内存≥64GB | BWA-MEM2 | DeepVariant CPU |
| `use_gpu=false` + 内存<64GB | BWA-MEM | DeepVariant CPU |

## License

MIT