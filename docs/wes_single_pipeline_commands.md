# 单样本全外显子测序(WES)分析流程命令

基于 Nextflow DSL2 的胚系变异检测流程。

## 目录

- [1. 质控过滤 (FASTP)](#1-质控过滤-fastp)
- [2. 序列比对 (BWA-MEM2)](#2-序列比对-bwa-mem2)
- [3. 标记重复 (GATK MarkDuplicates)](#3-标记重复-gatk-markduplicates)
- [4. 性别检查 (SEX_CHECK)](#4-性别检查-sex_check)
- [5. 覆盖度统计 (XAMDST)](#5-覆盖度统计-xamdst)
- [6. SNP/INDEL 变异检测 (DeepVariant)](#6-snpindel-变异检测-deepvariant)
- [7. 样本指纹验证 (SAMPLE_FINGERPRINT)](#7-样本指纹验证-sample_fingerprint)
- [8. 线粒体变异检测 (GATK Mutect2)](#8-线粒体变异检测-gatk-mutect2)
- [9. 单倍型定相 (WhatsHap)](#9-单倍型定相-whatshap)
- [10. 变异注释 (VEP)](#10-变异注释-vep)
- [11. 遗传模式注释 (GenMod)](#11-遗传模式注释-genmod)
- [12. STR 检测 (ExpansionHunter)](#12-str-检测-expansionhunter)
- [13. STR 注释 (Stranger)](#13-str-注释-stranger)
- [14. CNV 检测 (CNVkit)](#14-cnv-检测-cnvkit)
- [15. CNV 导出 VCF (CNVkit)](#15-cnv-导出-vcf-cnvkit)
- [16. 基因型分析 (PLINK2)](#16-基因型分析-plink2)
- [17. ROH 检测 (BCFtools)](#17-roh-检测-bcftools)
- [18. VCF 统计 (BCFtools)](#18-vcf-统计-bcftools)
- [19. 变异过滤 (Slivar)](#19-变异过滤-slivar)
- [20. 复合杂合子检测 (Slivar)](#20-复合杂合子检测-slivar)
- [21. TSV 导出 (Slivar)](#21-tsv-导出-slivar)
- [22. CNV/SV 注释 (SVDB)](#22-cnvsv-注释-svdb)

---

## 1. 质控过滤 (FASTP)

```bash
fastp \
    --in1 sample_R1.fq.gz \
    --in2 sample_R2.fq.gz \
    --out1 sample_R1.fq.gz \
    --out2 sample_R2.fq.gz \
    --json sample.fastp.json \
    --html sample.fastp.html \
    --thread ${threads} \
    --detect_adapter_for_pe
```

**说明**：
- 对 PE 测序数据进行质量控制和过滤
- 自动检测并去除 adapter
- 最多支持 16 线程

---

## 2. 序列比对 (BWA-MEM2)

### BWA-MEM2 (高内存版本)

```bash
bwa-mem2 mem \
    -t ${task.cpus} \
    -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
    ${fasta} \
    ${reads[0]} \
    ${reads[1]} | \
    samtools sort -@ ${task.cpus} --reference ${fasta} -O cram -o sample.sorted.cram -

samtools index -@ ${task.cpus} sample.sorted.cram
```

### BWA-MEM (低内存版本)

```bash
bwa mem \
    -t ${task.cpus} \
    -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
    ${fasta} \
    ${reads[0]} \
    ${reads[1]} | \
    samtools sort -@ ${task.cpus} --reference ${fasta} -O cram -o sample.sorted.cram -

samtools index -@ ${task.cpus} sample.sorted.cram
```

**说明**：
- BWA-MEM2 需要 >= 64GB 内存，否则自动降级为 BWA-MEM
- 添加 Read Group 信息
- 直接输出 CRAM 格式以节省存储空间

---

## 3. 标记重复 (GATK MarkDuplicates)

```bash
gatk MarkDuplicates \
    --INPUT sample.sorted.cram \
    --OUTPUT sample.md.cram \
    --METRICS_FILE sample.metrics.txt \
    --CREATE_INDEX false \
    --REFERENCE_SEQUENCE ${fasta}

# GATK 4.6.x 输出 CRAM 时不自动生成索引，使用 samtools
samtools index -@ ${task.cpus} sample.md.cram
```

**说明**：
- 标记 PCR 重复
- CRAM 格式输出时关闭 GATK 的 index 生成，改为 samtools 生成

---

## 4. 性别检查 (SEX_CHECK)

```bash
# 基于 SRY 基因覆盖度推断性别
# -F 2052: 排除未比对的 reads (0x4) 和重复 reads (0x400)
samtools view -c -F 2052 sample.cram "Y:2786355-2788241"
```

**说明**：
- 统计 SRY 基因区域的 reads 数（无需参考基因组）
- 男性样本 SRY 区域应有 reads（>= 10 判定为 male）
- 女性样本 SRY 区域几乎无 reads
- 支持 GRCh37 (Y:2654396-2656292) 和 GRCh38 (Y:2786355-2788241)

---

## 5. 覆盖度统计 (XAMDST)

```bash
xamdst \
    --target ${target_bed} \
    --bam ${alignment} \
    --outdir ./coverage

# 或使用 bedtools
bedtools coverage \
    -a ${target_bed} \
    -b ${alignment} \
    -d > sample.coverage.txt
```

**说明**：
- 计算目标区域的覆盖度统计
- 生成覆盖度分布报告

---

## 6. SNP/INDEL 变异检测 (DeepVariant)

```bash
/opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref=${fasta} \
    --reads=${alignment} \
    --output_vcf=sample.vcf.gz \
    --num_shards=${task.cpus} \
    --intermediate_results_dir=./tmp
```

**说明**：
- 基于深度学习的变异检测
- 使用 WES 模型
- GPU 加速效果更好

---

## 7. 样本指纹验证 (SAMPLE_FINGERPRINT)

```bash
bcftools isec \
    --prefix sample_fingerprint \
    -n=1 \
    -C \
    sample.vcf.gz \
    snp_panel.vcf.gz

# 计算匹配率
python check_fingerprint.py \
    sample_fingerprint/0000.vcf \
    snp_panel.vcf.gz
```

**说明**：
- 验证样本身份
- 检测样本混淆或污染

---

## 8. 线粒体变异检测 (GATK Mutect2)

```bash
gatk Mutect2 \
    --reference ${fasta} \
    --input ${alignment} \
    --output sample.mt.vcf.gz \
    --mitochondria-mode \
    --max-reads-per-alignment-start 0 \
    --max-mnp-distance 0 \
    -L MT,chrM,chrMT,M
```

**说明**：
- 线粒体模式自动优化参数
- 支持低频异质性变异检测
- 仅分析线粒体基因组区域

---

## 9. 单倍型定相 (WhatsHap)

```bash
whatshap phase \
    --indels \
    --reference ${fasta} \
    --output sample.phase.vcf \
    ${vcf} ${bam}

bgzip -f sample.phase.vcf
tabix -f -p vcf sample.phase.vcf.gz
```

**说明**：
- 基于 reads 比对信息进行单倍型分相
- 支持 indels phasing
- 输出带 HP 标签的 phased VCF

---

## 10. 变异注释 (VEP)

```bash
vep \
    --input_file ${vcf} \
    --output_file sample.vep.vcf.gz \
    --format vcf \
    --vcf \
    --compress_output bgzip \
    --offline \
    --cache \
    --dir_cache ${cache} \
    --dir_plugins ${plugins} \
    --refseq \
    --fasta ${fasta} \
    --assembly GRCh38 \
    --force_overwrite \
    --fork ${task.cpus} \
    --stats_file sample.vep_summary.html \
    --pick \
    --pick_order mane_plus_clinical,mane_select,refseq,canonical,ensembl \
    --max_af \
    --hgvs \
    --symbol \
    --uniprot \
    --domains \
    --pubmed

tabix -p vcf sample.vep.vcf.gz
```

**说明**：
- Ensembl VEP 变异功能注释
- 支持 RefSeq 转录本
- 自定义注释：ClinVar, InterVar, cytoBand, dbNSFP, SpliceAI
- pick 模式选择一个最佳转录本

---

## 11. 遗传模式注释 (GenMod)

```bash
genmod models \
    ${vcf} \
    --family_file ${ped_file} \
    --split_variants \
    --phased \
    --strict | bgzip -c > sample.genmod_models.vcf.gz

tabix -f -p vcf sample.genmod_models.vcf.gz
```

**说明**：
- 注释遗传继承模式 (AD, AR, X-linked 等)
- 复合杂合子检测
- CADD 分数注释
- 家系分析支持

---

## 12. STR 检测 (ExpansionHunter)

```bash
ExpansionHunter \
    --reads ${alignment} \
    --reference ${fasta} \
    --variant-catalog ${catalog.json} \
    --output-prefix sample \
    --sex ${sex} \
    --threads ${task.cpus}

bgzip -c sample.vcf > sample.vcf.gz
tabix -p vcf sample.vcf.gz
```

**说明**：
- 短串联重复序列 (STR) 扩展检测
- 检测与遗传病相关的 STR 扩展（如亨廷顿病、脆性X综合征）

---

## 13. STR 注释 (Stranger)

```bash
stranger \
    ${vcf} \
    --repeats-file ${repeats.json} | bgzip -c > sample.stranger.vcf.gz

tabix -f -p vcf sample.stranger.vcf.gz
```

**说明**：
- 注释 STR 重复序列大小
- 添加 STR_STATUS 字段 (normal/pre_mutation/full_mutation)

---

## 14. CNV 检测 (CNVkit)

```bash
# Step 1: 计算覆盖度
cnvkit.py coverage \
    ${alignment} \
    ${targets} \
    --fasta ${fasta} \
    --processes ${task.cpus} \
    --output sample.targetcoverage.cnn

cnvkit.py coverage \
    ${alignment} \
    ${antitargets} \
    --fasta ${fasta} \
    --processes ${task.cpus} \
    --output sample.antitargetcoverage.cnn

# Step 2: 修正覆盖度偏差
cnvkit.py fix \
    sample.targetcoverage.cnn \
    sample.antitargetcoverage.cnn \
    ${reference} \
    --output sample.cnr

# Step 3: CBS 分段
cnvkit.py segment \
    sample.cnr \
    --processes ${task.cpus} \
    --output sample.cns

# Step 4: CNV 调用
cnvkit.py call \
    sample.cns \
    --male-reference \
    --output sample.call.cns
```

**说明**：
- 基于测序深度的拷贝数变异检测
- coverage → fix → segment → call 四步流程

---

## 15. CNV 导出 VCF (CNVkit)

```bash
cnvkit.py export vcf \
    sample.call.cns \
    --male-reference \
    --output sample.cnv.vcf

bgzip -c sample.cnv.vcf > sample.cnv.vcf.gz
tabix -p vcf sample.cnv.vcf.gz
```

**说明**：
- 将 CNV 检测结果导出为 VCF 格式

---

## 16. 基因型分析 (PLINK2)

```bash
plink2 \
    --vcf sample.phase.vcf.gz \
    --make-pgen \
    --out sample

# 或输出 PLINK1 格式
plink2 \
    --vcf sample.phase.vcf.gz \
    --make-bed \
    --out sample
```

**说明**：
- VCF 转换为 PLINK2 格式
- 支持 PCA、亲缘关系等分析

---

## 17. ROH 检测 (BCFtools)

```bash
bcftools roh \
    --threads ${task.cpus} \
    --AF-dflt 0.001 \
    --output sample.roh \
    ${vcf}
```

**说明**：
- 检测基因组中的纯合子区域 (Runs of Homozygosity)
- 用于近亲繁殖分析、家系研究

---

## 18. VCF 统计 (BCFtools)

```bash
bcftools stats \
    --threads ${task.cpus} \
    --reference ${fasta} \
    ${vcf} > sample.vcfstats.txt
```

**说明**：
- 生成 VCF 文件的统计信息
- 包含变异数量、类型分布等

---

## 19. 变异过滤 (Slivar)

```bash
slivar expr \
    --vcf ${vcf} \
    --ped ${ped_file} \
    --info '${info_expr}' \
    --sample-expr '${sample_expr}' \
    --out-vcf sample.slivar.vcf.gz
```

**说明**：
- 使用 JavaScript 表达式进行 VCF 变异过滤
- 支持复杂过滤条件

---

## 20. 复合杂合子检测 (Slivar)

```bash
slivar compound-hets \
    --vcf ${vcf} \
    --ped ${ped_file} \
    --out-vcf sample.compound.vcf.gz
```

**说明**：
- 检测复合杂合子变异
- 两个不同位点的杂合变异

---

## 21. TSV 导出 (Slivar)

```bash
slivar tsv \
    --vcf ${vcf} \
    --ped ${ped_file} \
    --output sample.slivar.tsv.gz
```

**说明**：
- 将 VCF 导出为 TSV 表格
- 便于后续分析

---

## 22. CNV/SV 注释 (SVDB)

```bash
svdb \
    --query \
    --bnd_distance 2500 \
    --overlap 0.8 \
    --in_occ in_occ \
    --in_frq in_frq \
    --out_occ out_occ \
    --out_frq out_frq \
    --db ${svdb.vcf} \
    --vcf sample.cnv.vcf.gz \
    --out sample.cnv.svdb.vcf.gz
```

**说明**：
- CNV/SV 变异注释
- 添加数据库频率和计数信息

---

## 流程图

```
FASTP (质控过滤)
    ↓
BWA_MEM2 (比对)
    ↓
GATK_MARKDUPLICATES (标记重复) → samtools index (CRAM索引)
    ↓
    ├─→ SEX_CHECK (性别检查)
    ├─→ XAMDST (覆盖度统计)
    ├─→ DEEPVARIANT (SNP/INDEL检测)
    │       ↓
    │   WHATSHAP_PHASE (单倍型定相)
    │       ↓
    │   VEP (变异注释)
    │       ↓
    │   GENMOD_MODELS (遗传模式注释)
    │       ↓
    │   SLIVAR_EXPR (变异过滤)
    │       ├─→ SLIVAR_COMPOUND_HETS (复合杂合子)
    │       └─→ SLIVAR_TSV (TSV导出)
    │
    ├─→ GATK_MUTECT2_MT (线粒体变异)
    │       ↓
    │   GENMOD_MODELS_MT (线粒体遗传模式)
    │
    ├─→ EXPANSIONHUNTER (STR检测)
    │       ↓
    │   STRANGER (STR注释)
    │       ↓
    │   GENMOD_MODELS_STR (STR遗传模式)
    │
    └─→ CNVKIT_CALL (CNV检测)
            ↓
        CNVKIT_EXPORT_VCF (导出VCF)
            ↓
        GENMOD_MODELS_CNV (CNV遗传模式)
            ↓
        SVDB_QUERY (CNV/SV注释)

BCFTOOLS_ROH (ROH检测)
BCFTOOLS_STATS (VCF统计)
PLINK2_VCF_TO_PLINK (基因型分析)
```

---

## 资源需求

| 模块 | CPU | 内存 | 说明 |
|------|-----|------|------|
| FASTP | 8 | 4 GB | 多线程有效，内存需求低 |
| BWA_MEM2 | 16 | 64 GB | 内存 >= 64GB 用 bwa-mem2，否则降级为 bwa |
| GATK_MARKDUPLICATES | 4 | 16 GB | - |
| DEEPVARIANT | 16 | 32 GB | GPU 可加速 |
| VEP | 4 | 8 GB | - |
| EXPANSIONHUNTER | 4 | 8 GB | - |
| CNVKIT_CALL | 4 | 8 GB | - |

---

## 容器镜像

| 工具 | 镜像 |
|------|------|
| FASTP/BWA/Samtools/xamdst | `ghcr.io/pzweuj/mapping:2025dec` |
| GATK (含 samtools) | `ghcr.io/pzweuj/gatk:4.6.2.0` |
| DeepVariant | `google/deepvariant:1.6.0` |
| VEP | `ensemblorg/ensembl-vep:release_112.0` |
| ExpansionHunter | `clinicalgenomics/expansionhunter:5.0.0` |
| GenMod/WhatsHap/PLINK2 | `ghcr.io/pzweuj/mapping:2025dec` |
