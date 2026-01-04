#!/usr/bin/env python3
"""
Aggregate QC metrics from multiple sources into a comprehensive QC report.

Aggregates data from:
- fastp: JSON quality control reports
- bamdst: Coverage statistics
- GATK MarkDuplicates: Duplicate and alignment metrics
- Sex check: Sample sex verification
- VEP/DeepVariant: Variant calling stats (if available)

Usage:
    python aggregate_qc.py --input results/ --output qc_summary.tsv
    python aggregate_qc.py --input results/ --output qc_summary.tsv --json-output qc_summary.json
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Optional


def parse_fastp_json(json_file: Path) -> dict:
    """Parse fastp JSON report and extract key metrics."""
    try:
        data = json.loads(json_file.read_text())
        summary = data.get('summary', {})

        return {
            'sample_id': json_file.name.replace('.fastp.json', ''),
            'fastp_total_reads': summary.get('total_reads'),
            'fastp_total_bases': summary.get('total_bases'),
            'fastp_q30_rate': summary.get('q30_rate'),
            'fastp_gc_content': summary.get('gc_content'),
            'fastp_read1_adapter_detected': data.get('read1_after_filtering', {}).get('adapter_detected_seq', ''),
            'fastp_read2_adapter_detected': data.get('read2_after_filtering', {}).get('adapter_detected_seq', ''),
            'fastp_pct_reads_filtered': round(
                (1 - summary.get('total_reads', 0) / data.get('before_filtering', {}).get('total_reads', 1)) * 100
                if data.get('before_filtering', {}).get('total_reads', 0) > 0 else 0, 2
            )
        }
    except Exception as e:
        print(f"Warning: Failed to parse {json_file}: {e}", file=sys.stderr)
        return {}


def parse_bamdst_coverage(coverage_file: Path) -> dict:
    """Parse bamdst coverage.report and extract key metrics."""
    metrics = {}
    try:
        content = coverage_file.read_text()
        lines = content.strip().split('\n')

        for line in lines:
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip().lower().replace(' ', '_')
                value = value.strip()

                # Convert numeric values
                if value.replace('.', '').replace(',', '').isdigit():
                    if '.' in value:
                        value = float(value)
                    else:
                        value = int(value)
                elif '%' in value:
                    value = float(value.replace('%', ''))
                elif value in ['Yes', 'yes', 'TRUE', 'True']:
                    value = True
                elif value in ['No', 'no', 'FALSE', 'False']:
                    value = False

                metrics[f'bamdst_{key}'] = value

        # Extract sample_id from filename
        sample_id = coverage_file.name.replace('.coverage.stat', '').replace('.stat', '')
        if not metrics.get('bamdst_sample_id'):
            metrics['bamdst_sample_id'] = sample_id

    except Exception as e:
        print(f"Warning: Failed to parse {coverage_file}: {e}", file=sys.stderr)

    return metrics


def parse_gatk_metrics(metrics_file: Path) -> dict:
    """Parse GATK MarkDuplicates metrics file."""
    metrics = {}
    try:
        content = metrics_file.read_text()
        lines = content.strip().split('\n')

        # Find the data section
        in_data = False
        for line in lines:
            if line.startswith('## METRICS'):
                in_data = True
                continue
            if in_data and line.startswith('##'):
                break
            if in_data and line.strip():
                # Header line
                if 'LIBRARY' in line.upper():
                    headers = line.split('\t')
                # Data line
                elif 'Unknown' in line or 'all' in line.lower():
                    values = line.split('\t')
                    for i, header in enumerate(headers):
                        header = header.lower().replace(' ', '_')
                        if header in ['library', 'sample', 'read_group']:
                            continue
                        try:
                            val = float(values[i])
                            metrics[f'gatk_{header}'] = val
                        except (ValueError, IndexError):
                            pass
                    break

    except Exception as e:
        print(f"Warning: Failed to parse {metrics_file}: {e}", file=sys.stderr)

    return metrics


def parse_sex_check(txt_file: Path) -> dict:
    """Parse sex check output file."""
    metrics = {}
    try:
        content = txt_file.read_text()
        for line in content.strip().split('\n'):
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip().lower().replace(' ', '_')
                value = value.strip()
                metrics[f'sex_{key}'] = value
    except Exception as e:
        print(f"Warning: Failed to parse {txt_file}: {e}", file=sys.stderr)

    return metrics


def parse_deepvariant_stats(vcf_file: Path) -> dict:
    """Parse DeepVariant VCF for variant counts."""
    metrics = {}
    try:
        content = vcf_file.read_text()
        lines = content.strip().split('\n')

        snps = 0
        indels = 0
        for line in lines:
            if not line.startswith('#'):
                if 'INDEL' in line:
                    indels += 1
                else:
                    snps += 1

        sample_id = vcf_file.name.split('.')[0]
        metrics['dv_total_variants'] = snps + indels
        metrics['dv_snps'] = snps
        metrics['dv_indels'] = indels
        metrics['dv_sample_id'] = sample_id

    except Exception as e:
        print(f"Warning: Failed to parse {vcf_file}: {e}", file=sys.stderr)

    return metrics


def aggregate_qc_results(input_dir: Path) -> list[dict]:
    """Aggregate all QC results from the input directory."""
    results = {}

    # Parse fastp JSON files
    for json_file in input_dir.rglob('*fastp.json'):
        sample_metrics = parse_fastp_json(json_file)
        sample_id = sample_metrics.get('sample_id')
        if sample_id:
            results[sample_id] = results.get(sample_id, {})
            results[sample_id].update(sample_metrics)

    # Parse bamdst coverage files
    for cov_file in input_dir.rglob('*coverage.stat'):
        sample_metrics = parse_bamdst_coverage(cov_file)
        sample_id = sample_metrics.get('bamdst_sample_id')
        if sample_id:
            results[sample_id] = results.get(sample_id, {})
            results[sample_id].update(sample_metrics)

    # Parse GATK metrics files
    for metrics_file in input_dir.rglob('*.metrics.txt'):
        sample_metrics = parse_gatk_metrics(metrics_file)
        # Extract sample_id from filename pattern
        name = metrics_file.name
        if '.md.' in name:
            sample_id = name.split('.md.')[0]
        elif '.metrics' in name:
            sample_id = name.replace('.metrics.txt', '')
        else:
            sample_id = None

        if sample_id:
            results[sample_id] = results.get(sample_id, {})
            results[sample_id].update(sample_metrics)

    # Parse sex check files
    for sex_file in input_dir.rglob('*.sex_check.txt'):
        sample_metrics = parse_sex_check(sex_file)
        sample_id = sample_metrics.get('sex_sample_id')
        if sample_id:
            results[sample_id] = results.get(sample_id, {})
            results[sample_id].update(sample_metrics)

    # Parse DeepVariant VCF files for variant counts
    for vcf_file in input_dir.rglob('*.deepvariant.vcf.gz'):
        sample_metrics = parse_deepvariant_stats(vcf_file)
        sample_id = sample_metrics.get('dv_sample_id')
        if sample_id:
            results[sample_id] = results.get(sample_id, {})
            results[sample_id].update(sample_metrics)

    return list(results.values())


def write_tsv(results: list[dict], output_file: Path):
    """Write results to TSV format."""
    if not results:
        print("Warning: No results to write", file=sys.stderr)
        return

    # Get all unique keys
    all_keys = set()
    for r in results:
        all_keys.update(r.keys())

    # Remove duplicate sample_id columns
    all_keys = [k for k in all_keys if k != 'sample_id']
    all_keys.insert(0, 'sample_id')

    with open(output_file, 'w') as f:
        f.write('\t'.join(all_keys) + '\n')
        for r in results:
            row = []
            for key in all_keys:
                val = r.get(key, '')
                if val is None:
                    val = ''
                row.append(str(val))
            f.write('\t'.join(row) + '\n')


def write_json(results: list[dict], output_file: Path):
    """Write results to JSON format."""
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate QC metrics from multiple sources into a comprehensive report'
    )
    parser.add_argument('-i', '--input', required=True, type=Path,
                        help='Input directory containing QC results')
    parser.add_argument('-o', '--output', required=True, type=Path,
                        help='Output TSV file')
    parser.add_argument('--json-output', type=Path,
                        help='Optional JSON output file')
    parser.add_argument('--sample-prefix', default='',
                        help='Prefix to remove from sample IDs')

    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: Input directory not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"Aggregating QC results from {args.input}...")

    results = aggregate_qc_results(args.input)
    print(f"Found {len(results)} samples with QC data")

    # Write TSV output
    write_tsv(results, args.output)
    print(f"Written TSV report to {args.output}")

    # Write JSON output if requested
    if args.json_output:
        write_json(results, args.json_output)
        print(f"Written JSON report to {args.json_output}")


if __name__ == '__main__':
    main()
