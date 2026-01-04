#!/usr/bin/env python3
"""
Convert sample information to PLINK2 family format.

PLINK2 family (pedigree) format requirements:
- FAM_ID: Family ID
- IID: Individual ID
- PAT: Paternal ID (0 if unknown)
- MAT: Maternal ID (0 if unknown)
- SEX: 1=Male, 2=Female, 0=Unknown
- PHENOTYPE: -9=Missing, 0=Missing, 1=Unaffected, 2=Affected

Usage:
    python make_pedigree.py --input samples.yaml --output pedigree.fam
    python make_pedigree.py --input samples.json --output pedigree.fam
"""

import argparse
import json
import yaml
import sys
from pathlib import Path


def parse_sample_file(input_file: Path) -> list[dict]:
    """Parse sample information from YAML or JSON file."""
    content = input_file.read_text()

    if input_file.suffix in ['.yaml', '.yml']:
        samples = yaml.safe_load(content)
    elif input_file.suffix == '.json':
        samples = json.loads(content)
    else:
        # Try to auto-detect format
        try:
            samples = yaml.safe_load(content)
        except yaml.YAMLError:
            try:
                samples = json.loads(content)
            except json.JSONDecodeError:
                raise ValueError(f"Unable to parse {input_file}. Please use YAML or JSON format.")

    # Normalize to list
    if isinstance(samples, dict):
        # Single sample with sample_id as key
        if 'sample_id' in samples:
            samples = [samples]
        else:
            # Multiple samples with sample_id as keys
            samples = [{'sample_id': k, **v} for k, v in samples.items()]
    elif not isinstance(samples, list):
        raise ValueError(f"Expected list or dict, got {type(samples).__name__}")

    return samples


def samples_to_pedigree(samples: list[dict]) -> list[dict]:
    """Convert sample info to PLINK2 pedigree format."""
    pedigree = []

    # Build lookup for parental relationships
    sample_ids = {s['sample_id'] for s in samples}

    for sample in samples:
        sample_id = sample.get('sample_id')
        if not sample_id:
            print(f"Warning: Skipping sample without sample_id", file=sys.stderr)
            continue

        # Parse family ID (default to sample_id if not provided)
        fam_id = sample.get('family_id', sample_id)

        # Parse parental IDs
        father_id = sample.get('father_id')
        mother_id = sample.get('mother_id')

        # Convert to PLINK format (0 if unknown)
        pat = father_id if father_id and father_id in sample_ids else '0'
        mat = mother_id if mother_id and mother_id in sample_ids else '0'

        # Parse sex
        sex = sample.get('sex', 'unknown')
        sex_map = {
            'male': '1',
            'm': '1',
            'M': '1',
            '1': '1',
            'female': '2',
            'f': '2',
            'F': '2',
            '2': '2',
            'unknown': '0',
            '0': '0'
        }
        sex_plink = sex_map.get(str(sex).lower(), '0')

        # Parse phenotype
        phenotype = sample.get('phenotype', 'unknown')
        phenotype_map = {
            'affected': '2',
            '2': '2',
            'unaffected': '1',
            '1': '1',
            'unknown': '-9',
            'missing': '-9',
            '-9': '-9'
        }
        pheno = phenotype_map.get(str(phenotype).lower(), '-9')

        pedigree.append({
            'fam_id': fam_id,
            'iid': sample_id,
            'pat': pat,
            'mat': mat,
            'sex': sex_plink,
            'phenotype': pheno
        })

    return pedigree


def write_pedigree(pedigree: list[dict], output_file: Path):
    """Write pedigree to PLINK2 FAM file."""
    with open(output_file, 'w') as f:
        f.write('#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n')
        for ind in pedigree:
            f.write(f"{ind['fam_id']}\t{ind['iid']}\t{ind['pat']}\t{ind['mat']}\t"
                    f"{ind['sex']}\t{ind['phenotype']}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Convert sample information to PLINK2 family format'
    )
    parser.add_argument('-i', '--input', required=True, type=Path,
                        help='Input sample file (YAML or JSON)')
    parser.add_argument('-o', '--output', required=True, type=Path,
                        help='Output FAM file')
    parser.add_argument('--sample-key', default='samples',
                        help='Key for samples list (for JSON with nested structure)')

    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Parse samples
    samples = parse_sample_file(args.input)

    # Check for nested structure
    if isinstance(samples, dict) and args.sample_key in samples:
        samples = samples[args.sample_key]

    print(f"Loaded {len(samples)} samples from {args.input}")

    # Convert to pedigree format
    pedigree = samples_to_pedigree(samples)
    print(f"Converted {len(pedigree)} individuals to PLINK2 format")

    # Write output
    write_pedigree(pedigree, args.output)
    print(f"Written pedigree to {args.output}")


if __name__ == '__main__':
    main()
