# Sample Configuration File

## Configuration (sample.json)

All optional parameters included for reference.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sample_id` | string | Yes | Sample identifier |
| `read1` | string | Yes | Path to R1 FASTQ |
| `read2` | string | Yes | Path to R2 FASTQ |
| `reference.fasta` | string | Yes | Reference genome FASTA |
| `outdir` | string | Yes | Output directory |

### Optional Parameters

| Parameter | Description |
|-----------|-------------|
| `target_bed` | Target regions BED (for CNV/coverage) |
| `antitarget_bed` | Anti-target regions BED (for CNVkit) |
| `cnv_reference` | CNV reference baseline `.cnn` |
| `variant_catalog` | STR variant catalog JSON |
| `stranger_repeats` | STR repeat definitions JSON |
| `snp_panel` | Sample identity SNP panel VCF |
| `vep_cache` | VEP cache directory |
| `vep_plugins` | VEP plugin paths (array) |
| `vep_clinvar` | ClinVar VCF for annotation |
| `vep_intervar` | InterVar VCF for annotation |
| `slivar_js` | Slivar JavaScript functions |
| `gnotate` | gnomAD gnotate file |
| `svdb_db` | SVDB database VCF |
| `ped_file` | PED pedigree file |

## Usage

```bash
# 假设所有数据在 /data 目录，参考基因组在 /reference
WORKDIR=/data
REFDIR=/reference
RESULTS=/results

docker run --rm -it \
    -v ${WORKDIR}:${WORKDIR} \
    -v ${REFDIR}:${REFDIR} \
    -v ${RESULTS}:${RESULTS} \
    -w ${WORKDIR} \
    nextflow/nextflow:latest \
    nextflow run /path/to/schema-germline/main.nf \
        --config /path/to/examples/sample.json \
        -profile docker
```

> **重要**：JSON 中的所有路径必须是容器内可访问的绝对路径。建议将数据、参考基因组和结果目录挂载到容器内相同的路径。

## Reference File Requirements

All reference files must be pre-indexed:
- FASTA: `.fai` index
- BED: no index needed
- VCF: `.tbi` index (use `tabix`)
- BAM/CRAM: `.bai`/`.crai` index (use `samtools index`)

---

## CNV Baseline Pipeline

构建 CNVkit 参考基线，用于后续样本的 CNV 检测。

### CSV Format (cnv_baseline.csv)

```csv
sample_id,bam,bai
normal001,/path/to/normal001.cram,/path/to/normal001.cram.crai
normal002,/path/to/normal002.cram,/path/to/normal002.cram.crai
...
```

### JSON Config (cnv_baseline.json)

```json
{
  "samples": "/path/to/cnv_baseline.csv",
  "reference": {
    "fasta": "/path/to/reference/hg38.fa"
  },
  "target_bed": "/path/to/targets.bed",
  "annotate": "/path/to/reference/refFlat.txt",
  "access_bed": "/path/to/access.hg38.bed",
  "outdir": "/path/to/cnv_baseline"
}
```

### Parameters

| Parameter | Required | Description |
|-----------|----------|-------------|
| `samples` | Yes | Path to CSV file with sample list |
| `reference.fasta` | Yes | Reference genome FASTA |
| `target_bed` | Yes | Target regions BED |
| `annotate` | No | refFlat.txt for gene annotation |
| `access_bed` | No | Accessible regions BED |
| `outdir` | Yes | Output directory |

### Usage

```bash
# 假设所有数据在 /data 目录，参考基因组在 /reference
WORKDIR=/data
REFDIR=/reference
RESULTS=/results

docker run --rm -it \
    -v ${WORKDIR}:${WORKDIR} \
    -v ${REFDIR}:${REFDIR} \
    -v ${RESULTS}:${RESULTS} \
    -w ${WORKDIR} \
    nextflow/nextflow:latest \
    nextflow run /path/to/schema-germline/main.nf \
        -entry CNV_BASELINE \
        --config /path/to/examples/cnv_baseline.json \
        -profile docker
```

### Output

```
{outdir}/
├── targets.bed              # Target regions
├── antitargets.bed          # Anti-target regions
└── reference.cnn            # CNV reference baseline
```
