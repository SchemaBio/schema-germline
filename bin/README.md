# Utility Scripts for Schema Germline Pipeline

This directory contains utility scripts for data processing and QC reporting.

## make_pedigree.py

Convert sample information to PLINK2 family format.

### Input Format

**YAML format:**
```yaml
- sample_id: Sample001
  family_id: FAM001
  father_id: Sample002
  mother_id: Sample003
  sex: male
  phenotype: affected

- sample_id: Sample002
  family_id: FAM001
  sex: male
  phenotype: unaffected
```

**JSON format:**
```json
[
  {
    "sample_id": "Sample001",
    "family_id": "FAM001",
    "father_id": "Sample002",
    "mother_id": "Sample003",
    "sex": "male",
    "phenotype": "affected"
  }
]
```

### Output Format (PLINK2 FAM)

```
#FID    IID     PAT     MAT     SEX     PHENO
FAM001  Sample001       Sample002       Sample003       1       2
FAM001  Sample002       0       0       1       1
```

### Usage

```bash
python make_pedigree.py --input samples.yaml --output pedigree.fam
python make_pedigree.py --input samples.json --output pedigree.fam
```

### Options

| Option | Description |
|--------|-------------|
| `-i, --input` | Input sample file (YAML or JSON) |
| `-o, --output` | Output FAM file |
| `--sample-key` | Key for samples list in JSON (default: samples) |

---

## aggregate_qc.py

Aggregate QC metrics from multiple sources into a comprehensive TSV report.

### Input

Recursively searches the input directory for:
- `*fastp.json` - FASTP quality control reports
- `*coverage.stat` - Bamdst coverage statistics
- `*.metrics.txt` - GATK MarkDuplicates metrics
- `*.sex_check.txt` - Sex verification results
- `*.deepvariant.vcf.gz` - DeepVariant VCF files

### Output

Generates a TSV file with columns for:
- **FASTP**: total_reads, q30_rate, gc_content, pct_reads_filtered
- **Bamdst**: target_covered, mean_depth, coverage_rate_*
- **GATK**: duplication_rate, mapped_reads, pct_properly_paired
- **Sex Check**: inferred_sex, declared_sex, status
- **DeepVariant**: total_variants, snps, indels

### Usage

```bash
# Aggregate all QC results
python aggregate_qc.py --input results/ --output qc_summary.tsv

# Also generate JSON output
python aggregate_qc.py --input results/ --output qc_summary.tsv --json-output qc_summary.json
```

### Options

| Option | Description |
|--------|-------------|
| `-i, --input` | Input directory containing QC results |
| `-o, --output` | Output TSV file |
| `--json-output` | Optional JSON output file |
| `--sample-prefix` | Prefix to remove from sample IDs |

---

## Integration with Nextflow

These scripts can be used in Nextflow workflows:

```groovy
process AGGREGATE_QC {
    tag "qc_summary"
    label 'process_low'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path qc_dir

    output:
    path "qc_summary.tsv"

    script:
    """
    aggregate_qc.py --input ${qc_dir} --output qc_summary.tsv
    """
}
```
