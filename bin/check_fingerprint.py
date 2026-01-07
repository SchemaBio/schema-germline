#!/usr/bin/env python3
"""
Check sample fingerprint concordance against a SNP panel.

This script compares sample genotypes with expected genotypes from a SNP panel
to verify sample identity and detect sample mixups or contamination.

Usage:
    python check_fingerprint.py isec_output.vcf snp_panel.vcf.gz -o fingerprint_report.tsv
    python check_fingerprint.py isec_output.vcf snp_panel.vcf.gz --plot

Input:
    - isec_output.vcf: VCF from bcftools isec (sample genotypes at panel SNPs)
    - snp_panel.vcf.gz: SNP panel VCF with expected genotypes

Output:
    - Concordance rate and detailed mismatch information
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Optional


def parse_gt(gt_string: str) -> tuple:
    """
    Parse VCF genotype string to allele tuple.

    Args:
        gt_string: Genotype string like "0/1", "1/1", "./.", "0|1"

    Returns:
        Tuple of alleles (e.g., (0, 1), (1, 1), (None, None))
    """
    if gt_string == './.' or gt_string == '.|.':
        return (None, None)

    # Handle phased (|) and unphased (/) separators
    sep = '|' if '|' in gt_string else '/'
    parts = gt_string.split(sep)

    if len(parts) != 2:
        return (None, None)

    try:
        return (int(parts[0]), int(parts[1]))
    except ValueError:
        return (None, None)


def is_heterozygous(gt: tuple) -> bool:
    """Check if genotype is heterozygous (different alleles)."""
    if gt[0] is None or gt[1] is None:
        return False
    return gt[0] != gt[1] and gt[0] != -1 and gt[1] != -1


def is_homozygous(gt: tuple) -> bool:
    """Check if genotype is homozygous (same non-ref alleles)."""
    if gt[0] is None or gt[1] is None:
        return False
    return gt[0] == gt[1] and gt[0] != -1


def is_matching(gt1: tuple, gt2: tuple) -> Optional[bool]:
    """
    Check if two genotypes match.

    Returns:
        - True: genotypes match
        - False: genotypes mismatch
        - None: one or both genotypes are missing (./.)
    """
    if gt1[0] is None or gt2[0] is None:
        return None

    # Handle ref/ref (0/0) matching
    if gt1 == (0, 0) and gt2 == (0, 0):
        return True
    if gt1 == (0, 0) or gt2 == (0, 0):
        return False

    # For non-ref genotypes, check if they share any allele
    return set(gt1) & set(gt2)


def parse_vcf_samples(vcf_file: Path) -> dict:
    """
    Parse VCF file and extract sample genotypes.

    Args:
        vcf_file: Path to VCF file (can be gzipped)

    Returns:
        Dictionary mapping CHROM:POS to sample genotype
    """
    genotypes = {}

    import gzip

    open_func = gzip.open if str(vcf_file).endswith('.gz') else open

    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            chrom = parts[0]
            pos = parts[1]
            key = f"{chrom}:{pos}"

            # Get sample genotype (first sample column)
            sample_gt = parts[9].split(':')[0]  # Handle GT only or GT:AD format
            genotypes[key] = parse_gt(sample_gt)

    return genotypes


def parse_vcf_with_info(vcf_file: Path) -> dict:
    """
    Parse VCF file and extract sample genotypes with full INFO.

    Args:
        vcf_file: Path to VCF file (can be gzipped)

    Returns:
        Dictionary with position -> {gt, alleles, info}
    """
    import gzip

    variants = {}
    open_func = gzip.open if str(vcf_file).endswith('.gz') else open

    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]
            info = parts[7]

            # Get sample genotype
            sample_gt = parts[9].split(':')[0]
            gt = parse_gt(sample_gt)

            key = f"{chrom}:{pos}"
            variants[key] = {
                'gt': gt,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'info': info
            }

    return variants


def calculate_concordance(sample_genotypes: dict, panel_genotypes: dict) -> dict:
    """
    Calculate concordance between sample and panel genotypes.

    Args:
        sample_genotypes: Dict of sample genotypes at panel positions
        panel_genotypes: Dict of expected panel genotypes

    Returns:
        Dictionary with concordance statistics
    """
    total_sites = 0
    concordant = 0
    discordant = 0
    missing_sample = 0
    missing_panel = 0
    mismatches = []

    for pos, panel_gt in panel_genotypes.items():
        total_sites += 1

        if pos not in sample_genotypes:
            missing_sample += 1
            continue

        sample_gt = sample_genotypes[pos]

        if panel_gt[0] is None:
            missing_panel += 1
            continue

        if sample_gt[0] is None:
            missing_sample += 1
            continue

        match = is_matching(sample_gt, panel_gt)

        if match is True:
            concordant += 1
        elif match is False:
            discordant += 1
            mismatches.append({
                'position': pos,
                'sample_gt': '/'.join(str(a) for a in sample_gt if a is not None),
                'panel_gt': '/'.join(str(a) for a in panel_gt if a is not None),
                'sample_het': is_heterozygous(sample_gt),
                'panel_het': is_heterozygous(panel_gt)
            })

    # Calculate rates
    comparable = concordant + discordant
    concordance_rate = (concordant / comparable * 100) if comparable > 0 else 0
    call_rate = ((total_sites - missing_sample) / total_sites * 100) if total_sites > 0 else 0

    return {
        'total_sites': total_sites,
        'concordant': concordant,
        'discordant': discordant,
        'missing_sample': missing_sample,
        'missing_panel': missing_panel,
        'concordance_rate': concordance_rate,
        'call_rate': call_rate,
        'mismatches': mismatches
    }


def check_fingerprint(isec_vcf: Path, panel_vcf: Path, output: Optional[Path] = None) -> dict:
    """
    Main function to check sample fingerprint concordance.

    Args:
        isec_vcf: VCF file from bcftools isec (sample vs panel intersection)
        panel_vcf: SNP panel VCF file
        output: Optional output file for JSON report

    Returns:
        Dictionary with concordance results
    """
    print(f"Parsing sample genotypes from: {isec_vcf}")
    sample_genotypes = parse_vcf_with_info(isec_vcf)
    print(f"  Found {len(sample_genotypes)} sites in sample VCF")

    print(f"Parsing SNP panel genotypes from: {panel_vcf}")
    panel_genotypes = parse_vcf_with_info(panel_vcf)
    print(f"  Found {len(panel_genotypes)} sites in panel VCF")

    # Filter panel to only sites present in sample VCF
    panel_in_isec = {pos: gt for pos, gt in panel_genotypes.items() if pos in sample_genotypes}
    print(f"  Sites in intersection: {len(panel_in_isec)}")

    results = calculate_concordance(sample_genotypes, panel_in_isec)

    # Build JSON output
    report = {
        "sample_fingerprint": {
            "total_sites": results['total_sites'],
            "concordant": results['concordant'],
            "discordant": results['discordant'],
            "missing_in_sample": results['missing_sample'],
            "missing_in_panel": results['missing_panel'],
            "concordance_rate": round(results['concordance_rate'], 2),
            "call_rate": round(results['call_rate'], 2),
            "status": "PASS" if results['concordance_rate'] >= 99 else ("WARNING" if results['concordance_rate'] >= 95 else "FAIL"),
            "mismatches": results['mismatches']
        }
    }

    # Print JSON to stdout
    print(json.dumps(report, indent=2))

    # Write JSON report
    if output:
        with open(output, 'w') as f:
            json.dump(report, f, indent=2)
        print(f"\nReport written to: {output}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Check sample fingerprint concordance against a SNP panel'
    )
    parser.add_argument('isec_vcf', type=Path,
                        help='VCF from bcftools isec (sample vs panel intersection)')
    parser.add_argument('panel_vcf', type=Path,
                        help='SNP panel VCF file with expected genotypes')
    parser.add_argument('-o', '--output', type=Path,
                        help='Output TSV file for detailed mismatch report')
    parser.add_argument('--sample-id', default='sample',
                        help='Sample ID for the report (default: sample)')

    args = parser.parse_args()

    if not args.isec_vcf.exists():
        print(f"Error: isec VCF not found: {args.isec_vcf}", file=sys.stderr)
        sys.exit(1)

    if not args.panel_vcf.exists():
        print(f"Error: Panel VCF not found: {args.panel_vcf}", file=sys.stderr)
        sys.exit(1)

    results = check_fingerprint(args.isec_vcf, args.panel_vcf, args.output)

    # Exit with appropriate code (for pipeline integration)
    if results['concordance_rate'] < 95:
        sys.exit(1)
    elif results['concordance_rate'] < 99:
        sys.exit(1)  # WARNING also fails for strict pipelines
    sys.exit(0)


if __name__ == '__main__':
    main()
