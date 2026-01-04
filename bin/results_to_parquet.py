#!/usr/bin/env python3
"""
Aggregate multiple result files into a single Parquet file.

This script collects all sample result files (QC metrics, variant calls,
coverage reports, etc.) and combines them into a structured Parquet file
with nested tables for each sample.

Structure:
- Each row represents one sample
- Sample-level fields are stored as columns
- Each sample's detailed results are stored as nested dictionaries

Usage:
    python results_to_parquet.py --input results/ --output samples.parquet
    python results_to_parquet.py --input results/ --output samples.parquet --sample-pattern "Sample_*"
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Optional


def parse_key_value_file(file_path: Path) -> dict[str, Any]:
    """Parse key:value format text file."""
    data = {}
    try:
        for line in file_path.read_text().strip().split('\n'):
            line = line.strip()
            if not line or ':' not in line:
                continue
            key, value = line.split(':', 1)
            key = key.strip().lower().replace(' ', '_').replace('-', '_')
            value = value.strip()

            # Convert value types
            if value.lower() in ['true', 'yes', 'y']:
                value = True
            elif value.lower() in ['false', 'no', 'n']:
                value = False
            elif value.lower() in ['null', 'none', 'na', 'n/a']:
                value = None
            elif value.replace('.', '').replace('-', '').isdigit():
                if '.' in value:
                    value = float(value)
                else:
                    value = int(value)
            elif value.endswith('%') and value[:-1].replace('.', '').isdigit():
                value = float(value[:-1]) / 100

            data[key] = value
    except Exception as e:
        print(f"Warning: Failed to parse {file_path}: {e}", file=sys.stderr)
    return data


def parse_delimited_file(file_path: Path, delimiter: str = '\t') -> list[dict]:
    """Parse tab/space/comma delimited file to list of dicts."""
    rows = []
    try:
        lines = file_path.read_text().strip().split('\n')
        if not lines:
            return rows

        # Parse header
        header = lines[0].strip().lstrip('#').split(delimiter)
        header = [h.strip().lower().replace(' ', '_').replace('-', '_') for h in header]

        # Parse data rows
        for line in lines[1:]:
            if line.startswith('#'):
                continue
            values = line.split(delimiter)
            if len(values) == len(header):
                row = {}
                for i, h in enumerate(header):
                    val = values[i].strip()
                    # Type conversion
                    if val.lower() in ['null', 'none', 'na']:
                        val = None
                    elif val.replace('.', '').replace('-', '').isdigit():
                        val = float(val) if '.' in val else int(val)
                    row[h] = val
                rows.append(row)
    except Exception as e:
        print(f"Warning: Failed to parse {file_path}: {e}", file=sys.stderr)
    return rows


def parse_json_file(file_path: Path) -> dict | list:
    """Parse JSON file."""
    try:
        return json.loads(file_path.read_text())
    except Exception as e:
        print(f"Warning: Failed to parse {file_path}: {e}", file=sys.stderr)
        return {}


def parse_vcf_summary(file_path: Path) -> dict:
    """Parse VCF file for variant summary."""
    summary = {
        'total_variants': 0,
        'snps': 0,
        'indels': 0,
        'mnps': 0,
        'insertions': 0,
        'deletions': 0,
        'het_variants': 0,
        'hom_variants': 0
    }
    try:
        for line in file_path.read_text().strip().split('\n'):
            if line.startswith('#'):
                continue
            summary['total_variants'] += 1
            if 'INDEL' in line:
                summary['indels'] += 1
            if 'VT=SNP' in line or ('DP=' in line and 'INDEL' not in line):
                summary['snps'] += 1
            if 'VT=MNP' in line:
                summary['mnps'] += 1
            if '0/1' in line or '1/0' in line:
                summary['het_variants'] += 1
            elif '1/1' in line or '2/2' in line:
                summary['hom_variants'] += 1
    except Exception as e:
        print(f"Warning: Failed to parse VCF {file_path}: {e}", file=sys.stderr)
    return summary


def parse_fastp_json(file_path: Path) -> dict:
    """Parse fastp JSON report."""
    try:
        data = json.loads(file_path.read_text())
        summary = data.get('summary', {})

        # Parse adapter sequences
        read1_adapter = data.get('read1_after_filtering', {}).get('adapter_detected_seq', '')
        read2_adapter = data.get('read2_after_filtering', {}).get('adapter_detected_seq', '')

        before = data.get('before_filtering', {})
        after = data.get('after_filtering', {})

        return {
            'total_reads_raw': before.get('total_reads'),
            'total_bases_raw': before.get('total_bases'),
            'total_reads_filtered': after.get('total_reads'),
            'q30_rate': summary.get('q30_rate'),
            'gc_content': summary.get('gc_content'),
            'adapter_detected': bool(read1_adapter or read2_adapter),
            'adapter_sequence': f"{read1_adapter};{read2_adapter}" if read1_adapter or read2_adapter else '',
            'pct_reads_retained': round(
                after.get('total_reads', 0) / before.get('total_reads', 1) * 100, 2
            ) if before.get('total_reads', 0) > 0 else 0
        }
    except Exception as e:
        print(f"Warning: Failed to parse fastp JSON {file_path}: {e}", file=sys.stderr)
        return {}


def extract_sample_id(filename: str, patterns: list[str] = None) -> Optional[str]:
    """Extract sample ID from filename using common patterns."""
    import re

    # Common patterns for sample naming
    sample_patterns = [
        r'^([A-Za-z0-9_-]+)[._]',      # SampleID_xxx or SampleID.xxx
        r'^([A-Za-z0-9_-]+)-[A-Za-z]+$',  # SampleID-toolname
        r'([A-Za-z0-9]+)\.[a-zA-Z]+$',    # SampleID.ext
    ]

    if patterns:
        sample_patterns = patterns + sample_patterns

    for pattern in sample_patterns:
        match = re.match(pattern, filename)
        if match:
            return match.group(1)

    # Return filename without extension as fallback
    return Path(filename).stem.split('.')[0]


def categorize_file(file_path: Path) -> str:
    """Categorize file type based on name and content."""
    name = file_path.name.lower()

    if 'fastp' in name and name.endswith('.json'):
        return 'fastp_json'
    elif 'coverage' in name and ('stat' in name or 'report' in name):
        return 'coverage'
    elif 'metrics' in name and 'markdup' in name:
        return 'markdup_metrics'
    elif 'sex' in name and 'check' in name:
        return 'sex_check'
    elif 'roh' in name and not name.endswith('.vcf'):
        return 'roh'
    elif 'vcfstats' in name or name.endswith('.vcfstats.txt'):
        return 'vcf_stats'
    elif name.endswith('.vcf.gz') or name.endswith('.vcf'):
        return 'vcf'
    elif 'distribution' in name:
        return 'depth_distribution'
    elif 'chromosome' in name and 'stat' in name:
        return 'chromosome_stat'
    elif 'region' in name and 'stat' in name:
        return 'region_stat'
    elif 'genelist' in name and 'stat' in name:
        return 'genelist_stat'
    elif name.endswith('.json'):
        return 'json'
    elif name.endswith('.yaml') or name.endswith('.yml'):
        return 'yaml'
    else:
        return 'other'


def process_directory(input_dir: Path, sample_pattern: str = None) -> dict[str, dict]:
    """Process directory and organize results by sample."""
    samples = {}

    for file_path in input_dir.rglob('*'):
        if not file_path.is_file():
            continue

        # Skip hidden files and directories
        if file_path.name.startswith('.'):
            continue

        # Extract sample ID
        sample_id = extract_sample_id(file_path.name)
        if not sample_id:
            continue

        # Filter by pattern if specified
        if sample_pattern and not sample_id.startswith(sample_pattern.replace('*', '')):
            continue

        # Initialize sample dict
        if sample_id not in samples:
            samples[sample_id] = {'_meta': {'sample_id': sample_id, 'files': []}}

        # Categorize and parse file
        category = categorize_file(file_path)
        file_info = {
            'path': str(file_path),
            'category': category,
            'filename': file_path.name
        }

        samples[sample_id]['_meta']['files'].append(file_info)

        # Parse file content based on category
        if category == 'fastp_json':
            samples[sample_id]['fastp'] = parse_fastp_json(file_path)
        elif category == 'coverage':
            cov_data = parse_key_value_file(file_path)
            samples[sample_id]['coverage'] = cov_data
        elif category == 'markdup_metrics':
            dup_data = parse_key_value_file(file_path)
            samples[sample_id]['markdup'] = dup_data
        elif category == 'sex_check':
            sex_data = parse_key_value_file(file_path)
            samples[sample_id]['sex_check'] = sex_data
        elif category == 'roh':
            roh_data = parse_key_value_file(file_path)
            samples[sample_id]['roh'] = roh_data
        elif category == 'vcf_stats':
            vcf_stats = parse_delimited_file(file_path)
            samples[sample_id]['vcf_stats'] = vcf_stats
        elif category == 'vcf':
            vcf_summary = parse_vcf_summary(file_path)
            samples[sample_id]['variant_summary'] = vcf_summary
        elif category == 'depth_distribution':
            dist_data = parse_delimited_file(file_path)
            samples[sample_id]['depth_distribution'] = dist_data
        elif category == 'chromosome_stat':
            chr_data = parse_delimited_file(file_path)
            samples[sample_id]['chromosome_coverage'] = chr_data
        elif category == 'region_stat':
            region_data = parse_delimited_file(file_path)
            samples[sample_id]['region_coverage'] = region_data
        elif category == 'genelist_stat':
            gene_data = parse_delimited_file(file_path)
            samples[sample_id]['gene_coverage'] = gene_data
        elif category == 'json':
            json_data = parse_json_file(file_path)
            if json_data:
                samples[sample_id].setdefault('json_data', {})[file_path.stem] = json_data
        elif category == 'other':
            # Try to parse as key-value
            kv_data = parse_key_value_file(file_path)
            if kv_data:
                samples[sample_id].setdefault('other', {})[file_path.stem] = kv_data
            else:
                # Try delimited
                delim_data = parse_delimited_file(file_path)
                if delim_data:
                    samples[sample_id].setdefault('other', {})[file_path.stem] = delim_data

    return samples


def samples_to_dataframe(samples: dict) -> dict:
    """Convert samples dict to flat and nested DataFrames."""
    import pandas as pd

    # Flat table: sample-level info
    flat_rows = []

    # Nested tables
    tables = {
        'samples': None,
        'fastp': [],
        'coverage': [],
        'markdup': [],
        'sex_check': [],
        'roh': [],
        'variant_summary': [],
        'chromosome_coverage': [],
        'region_coverage': [],
        'gene_coverage': [],
        'depth_distribution': []
    }

    for sample_id, data in samples.items():
        row = {'sample_id': sample_id}
        nested_tables = {}

        for key, value in data.items():
            if key == '_meta':
                row.update({f'meta_{k}': v for k, v in value.items() if k != 'files'})
            elif isinstance(value, dict) and not value:
                continue
            elif isinstance(value, dict):
                # Check if it's a flat dict (sample-level metrics)
                if all(not isinstance(v, (dict, list)) for v in value.values()):
                    row[key] = json.dumps(value) if value else None
                else:
                    nested_tables[key] = value
            elif isinstance(value, list):
                nested_tables[key] = value
            else:
                row[key] = value

        flat_rows.append(row)

        # Populate nested tables
        if 'fastp' in nested_tables:
            tables['fastp'].append({'sample_id': sample_id, **nested_tables['fastp']})
        if 'coverage' in nested_tables:
            tables['coverage'].append({'sample_id': sample_id, **nested_tables['coverage']})
        if 'markdup' in nested_tables:
            tables['markdup'].append({'sample_id': sample_id, **nested_tables['markdup']})
        if 'sex_check' in nested_tables:
            tables['sex_check'].append({'sample_id': sample_id, **nested_tables['sex_check']})
        if 'roh' in nested_tables:
            tables['roh'].append({'sample_id': sample_id, **nested_tables['roh']})
        if 'variant_summary' in nested_tables:
            tables['variant_summary'].append({'sample_id': sample_id, **nested_tables['variant_summary']})

        # Lists (chromosome, region, gene coverage)
        for key in ['chromosome_coverage', 'region_coverage', 'gene_coverage', 'depth_distribution']:
            if key in nested_tables and isinstance(nested_tables[key], list):
                for item in nested_tables[key]:
                    item['sample_id'] = sample_id
                    tables[key].append(item)

    # Create DataFrames
    df_samples = pd.DataFrame(flat_rows)

    for key in tables:
        if tables[key]:
            tables[key] = pd.DataFrame(tables[key])

    return df_samples, tables


def write_parquet(
    df_samples: "pd.DataFrame",
    tables: dict,
    output_file: Path,
    nested: bool = True
):
    """Write DataFrames to Parquet file.

    Args:
        df_samples: Main samples DataFrame
        tables: Dictionary of additional tables
        output_file: Output Parquet file path
        nested: If True, store nested data as JSON columns; if False, use separate tables
    """
    import pandas as pd

    if nested:
        # Convert dict/list columns to JSON strings for storage
        for col in df_samples.columns:
            if df_samples[col].dtype == 'object':
                # Keep as-is for now, pandas will handle
                pass

        # Write main table with all data
        df_samples.to_parquet(output_file, index=False, engine='pyarrow')

        # Write separate parquet files for nested tables
        base_name = output_file.stem
        for key, df in tables.items():
            if df is not None and not df.empty:
                table_file = output_file.parent / f"{base_name}_{key}.parquet"
                df.to_parquet(table_file, index=False, engine='pyarrow')
                print(f"Written {table_file}")
    else:
        # Single file with all data as nested columns
        df_samples.to_parquet(output_file, index=False, engine='pyarrow')

    print(f"Written {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate multiple result files into a single Parquet file'
    )
    parser.add_argument('-i', '--input', required=True, type=Path,
                        help='Input directory containing result files')
    parser.add_argument('-o', '--output', required=True, type=Path,
                        help='Output Parquet file')
    parser.add_argument('--sample-pattern', '-p',
                        help='Pattern to filter sample IDs (e.g., "Sample_*")')
    parser.add_argument('--single-table', action='store_true',
                        help='Store all data in single table (no separate nested tables)')
    parser.add_argument('--engine', default='pyarrow', choices=['pyarrow', 'fastparquet'],
                        help='Parquet engine to use')

    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: Input directory not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing results from {args.input}...")

    # Process directory
    samples = process_directory(args.input, args.sample_pattern)
    print(f"Found {len(samples)} samples")

    if not samples:
        print("Warning: No samples found", file=sys.stderr)
        sys.exit(0)

    # Convert to DataFrames
    df_samples, tables = samples_to_dataframe(samples)

    # Show summary
    print(f"\nSamples table: {len(df_samples)} rows, {len(df_samples.columns)} columns")
    for key, df in tables.items():
        if df is not None and not df.empty:
            print(f"  {key}: {len(df)} rows")

    # Write output
    write_parquet(df_samples, tables, args.output, nested=not args.single_table)

    print(f"\nOutput written to: {args.output}")


if __name__ == '__main__':
    main()
