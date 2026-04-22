#!/usr/bin/env python3
"""
Sample Identity Fingerprint Generator

Generate a unique identity fingerprint from SNP genotypes using bcftools.
This script uses the Pengelly SNP panel (24 SNPs) for sample tracking.
"""

import argparse
import subprocess
import tempfile
import os
import sys
from pathlib import Path


def load_snp_positions(snp_file: str, assembly: str = "grch38") -> list:
    """
    Load SNP positions from Pengelly SNP file.

    Args:
        snp_file: Path to pengelly_snp.txt
        assembly: "grch37" or "grch38"

    Returns:
        List of tuples: (chromosome, position, rsid, gene)
    """
    col_map = {"grch37": 2, "grch38": 3}
    pos_col = col_map.get(assembly.lower(), 3)

    positions = []
    with open(snp_file, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                chrom = parts[1]
                pos = parts[pos_col]
                rsid = parts[4]
                gene = parts[5]
                positions.append((chrom, pos, rsid, gene))

    return positions


def create_position_file(positions: list, output_file: str):
    """
    Create a position file for bcftools mpileup.

    Format: chr:position (1-based)
    """
    with open(output_file, 'w') as f:
        for chrom, pos, rsid, gene in positions:
            f.write(f"{chrom}:{pos}\n")


def run_bcftools_mpileup(bam: str, fasta: str, positions_file: str,
                          output_vcf: str, threads: int = 1,
                          max_dp: int = 200) -> str:
    """
    Run bcftools mpileup and call to generate VCF.

    Args:
        bam: Input BAM file
        fasta: Reference FASTA file
        positions_file: File with positions (chr:pos format)
        output_vcf: Output VCF file path
        threads: Number of threads
        min_dp: Minimum depth filter

    Returns:
        Path to output VCF
    """
    # Run mpileup with AD (allele depth) annotation
    mpileup_cmd = [
        "bcftools", "mpileup",
        "-f", fasta,
        "-l", positions_file,
        "-a", "AD,DP",
        "-d", str(max_dp),  # Limit max depth for efficiency
        "--threads", str(threads),
        bam
    ]

    # Run call with genotype likelihoods
    call_cmd = [
        "bcftools", "call",
        "-mv",
        "-Ov",
        "-o", output_vcf
    ]

    # Pipe mpileup to call
    mpileup_proc = subprocess.Popen(
        mpileup_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    call_proc = subprocess.Popen(
        call_cmd,
        stdin=mpileup_proc.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    mpileup_proc.stdout.close()

    stdout, stderr = call_proc.communicate()

    if call_proc.returncode != 0:
        print(f"Error running bcftools: {stderr.decode()}", file=sys.stderr)
        sys.exit(1)

    return output_vcf


def extract_genotypes(vcf_file: str, positions: list,
                       min_dp: int = 10, min_af: float = 0.15,
                       max_af_het: float = 0.85) -> dict:
    """
    Extract genotypes from VCF file with allele frequency-based determination.

    Args:
        vcf_file: Input VCF
        positions: List of (chrom, pos, rsid, gene)
        min_dp: Minimum depth threshold (default: 10)
        min_af: Minimum allele frequency to count as ALT (default: 0.15)
        max_af_het: Maximum AF for heterozygous (above this is homozygous ALT)

    Returns:
        Dictionary: rsid -> genotype info dict
    """
    # Create lookup by position
    pos_lookup = {(chrom, pos): (rsid, gene) for chrom, pos, rsid, gene in positions}

    genotypes = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alt = parts[4]
            info = parts[7]
            sample = parts[9] if len(parts) > 9 else ""

            key = (chrom, pos)
            if key in pos_lookup:
                rsid, gene = pos_lookup[key]

                # Parse FORMAT fields
                format_fields = parts[8].split(':')
                sample_values = sample.split(':')
                field_dict = dict(zip(format_fields, sample_values))

                # Get DP (total depth)
                dp = int(field_dict.get('DP', 0))

                # Get AD (allele depths: ref_depth,alt_depth1,alt_depth2,...)
                ad_str = field_dict.get('AD', '0')
                ad_values = [int(x) for x in ad_str.split(',')]

                # Calculate allele frequencies
                ref_depth = ad_values[0] if len(ad_values) > 0 else 0
                alt_depths = ad_values[1:] if len(ad_values) > 1 else []

                # Determine genotype based on frequencies
                genotype, status = determine_genotype(
                    ref, alt.split(','), ref_depth, alt_depths, dp,
                    min_dp, min_af, max_af_het
                )

                # Calculate allele frequencies for output
                af_values = [d/dp if dp > 0 else 0 for d in ad_values]

                genotypes[rsid] = {
                    'gene': gene,
                    'chrom': chrom,
                    'pos': pos,
                    'genotype': genotype,
                    'ref': ref,
                    'alt': alt,
                    'dp': dp,
                    'ad': ad_str,
                    'af': ','.join(f"{af:.2f}" for af in af_values),
                    'status': status  # 'confident', 'low_depth', 'ambiguous'
                }

    return genotypes


def determine_genotype(ref: str, alts: list, ref_depth: int,
                        alt_depths: list, dp: int,
                        min_dp: int, min_af: float, max_af_het: float) -> tuple:
    """
    Determine genotype based on allele frequencies.

    Args:
        ref: Reference allele
        alts: List of alternate alleles
        ref_depth: Depth of reference allele
        alt_depths: Depths of alternate alleles
        dp: Total depth
        min_dp: Minimum depth threshold
        min_af: Minimum AF to count as ALT
        max_af_het: Max AF for heterozygous

    Returns:
        (genotype_string, status)
    """
    if dp < min_dp:
        return "./.", "low_depth"

    alleles = [ref] + alts
    depths = [ref_depth] + alt_depths

    # Calculate frequencies
    frequencies = [(d/dp, a, d) for d, a in zip(depths, alleles)]

    # Sort by frequency (descending)
    frequencies.sort(key=lambda x: x[0], reverse=True)

    # Get top alleles above threshold
    significant_alleles = []
    for af, allele, depth in frequencies:
        if af >= min_af and depth > 0:
            significant_alleles.append((af, allele))

    # Determine genotype
    if len(significant_alleles) == 0:
        # No significant alleles - likely all reference with sequencing noise
        return ref + ref, "confident"

    elif len(significant_alleles) == 1:
        af, allele = significant_alleles[0]
        if af >= max_af_het:
            # Homozygous ALT
            return allele + allele, "confident"
        else:
            # Heterozygous with reference
            if allele != ref:
                # Sort for consistent output
                gt = ''.join(sorted([ref, allele]))
                return gt, "confident"
            else:
                return ref + ref, "confident"

    elif len(significant_alleles) == 2:
        # Two significant alleles - could be heterozygous or ambiguous
        af1, allele1 = significant_alleles[0]
        af2, allele2 = significant_alleles[1]

        # Check if it's a clean heterozygous (~50%/50%)
        ratio = min(af1, af2) / max(af1, af2)
        if ratio > 0.5:  # Roughly equal proportions
            gt = ''.join(sorted([allele1, allele2]))
            return gt, "confident"
        else:
            gt = ''.join(sorted([allele1, allele2]))
            return gt, "ambiguous"

    else:
        # Multiple alleles - ambiguous
        top2 = sorted([significant_alleles[0][1], significant_alleles[1][1]])
        return ''.join(top2), "ambiguous"


def parse_depth(sample: str, info: str) -> int:
    """Extract depth from VCF sample or info field."""
    parts = sample.split(':')

    # Try to find DP in sample fields
    # Common format: GT:AD:DP:GQ:PL
    if len(parts) >= 3:
        try:
            # DP is usually 3rd field after GT and AD
            dp = int(parts[2])
            return dp
        except ValueError:
            pass

    # Try to extract from INFO
    for item in info.split(';'):
        if item.startswith('DP='):
            try:
                return int(item.split('=')[1])
            except ValueError:
                pass

    return 0


def generate_fingerprint(genotypes: dict) -> str:
    """
    Generate a fingerprint string from genotypes.

    Args:
        genotypes: Dictionary of rsid -> genotype info

    Returns:
        Fingerprint string
    """
    # Sort by rsid for consistent ordering
    sorted_rsids = sorted(genotypes.keys(), key=lambda x: int(x.replace('rs', '')))

    fingerprint_parts = []
    for rsid in sorted_rsids:
        gt = genotypes[rsid]['genotype']
        fingerprint_parts.append(f"{rsid}={gt}")

    return ";".join(fingerprint_parts)


def generate_hash_fingerprint(genotypes: dict) -> str:
    """
    Generate a hash-based fingerprint.

    More compact representation suitable for database storage.
    """
    import hashlib

    fingerprint = generate_fingerprint(genotypes)
    return hashlib.sha256(fingerprint.encode()).hexdigest()[:16]


def write_output(genotypes: dict, output_file: str, format: str = "tsv"):
    """
    Write genotype output in specified format.

    Args:
        genotypes: Dictionary of genotypes
        output_file: Output file path
        format: "tsv", "json", or "fingerprint"
    """
    sorted_rsids = sorted(genotypes.keys(), key=lambda x: int(x.replace('rs', '')))

    with open(output_file, 'w') as f:
        if format == "tsv":
            f.write("rsid\tgene\tchrom\tpos\tgenotype\tref\talt\tdp\tad\taf\tstatus\n")
            for rsid in sorted_rsids:
                g = genotypes[rsid]
                f.write(f"{rsid}\t{g['gene']}\t{g['chrom']}\t{g['pos']}\t"
                        f"{g['genotype']}\t{g['ref']}\t{g['alt']}\t{g['dp']}\t"
                        f"{g['ad']}\t{g['af']}\t{g['status']}\n")

        elif format == "json":
            import json
            output = {
                "fingerprint": generate_fingerprint(genotypes),
                "fingerprint_hash": generate_hash_fingerprint(genotypes),
                "genotypes": {rsid: genotypes[rsid] for rsid in sorted_rsids}
            }
            json.dump(output, f, indent=2)

        elif format == "fingerprint":
            f.write(f"# Fingerprint: {generate_fingerprint(genotypes)}\n")
            f.write(f"# Hash: {generate_hash_fingerprint(genotypes)}\n")
            f.write("# Genotypes:\n")
            for rsid in sorted_rsids:
                g = genotypes[rsid]
                f.write(f"{rsid}\t{g['genotype']}\t{g['dp']}\t{g['status']}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate sample identity fingerprint from SNP genotypes"
    )
    parser.add_argument("-b", "--bam", required=True,
                        help="Input BAM file")
    parser.add_argument("-f", "--fasta", required=True,
                        help="Reference FASTA file")
    parser.add_argument("-s", "--snp-file", required=True,
                        help="Pengelly SNP positions file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file")
    parser.add_argument("-a", "--assembly", default="grch38",
                        choices=["grch37", "grch38"],
                        help="Reference assembly version (default: grch38)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (default: 1)")
    parser.add_argument("--format", default="tsv",
                        choices=["tsv", "json", "fingerprint"],
                        help="Output format (default: tsv)")
    parser.add_argument("--min-dp", type=int, default=10,
                        help="Minimum depth threshold (default: 10)")
    parser.add_argument("--max-dp", type=int, default=200,
                        help="Maximum depth to analyze per position (default: 200)")
    parser.add_argument("--min-af", type=float, default=0.15,
                        help="Minimum allele frequency to count as ALT (default: 0.15)")
    parser.add_argument("--max-af-het", type=float, default=0.85,
                        help="Maximum AF for heterozygous; above this is homozygous ALT (default: 0.85)")
    parser.add_argument("--keep-vcf", action="store_true",
                        help="Keep intermediate VCF file")

    args = parser.parse_args()

    # Load SNP positions
    print(f"Loading SNP positions from {args.snp_file}...")
    positions = load_snp_positions(args.snp_file, args.assembly)
    print(f"Loaded {len(positions)} SNP positions")

    # Create temporary position file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
        create_position_file(positions, tmp.name)
        positions_file = tmp.name

    try:
        # Run bcftools
        vcf_file = args.output + ".vcf"
        print(f"Running bcftools mpileup on {args.bam}...")
        run_bcftools_mpileup(args.bam, args.fasta, positions_file, vcf_file,
                              args.threads, args.max_dp)

        # Extract genotypes
        print("Extracting genotypes...")
        genotypes = extract_genotypes(
            vcf_file, positions,
            min_dp=args.min_dp,
            min_af=args.min_af,
            max_af_het=args.max_af_het
        )
        print(f"Extracted genotypes for {len(genotypes)} positions")

        # Count status
        status_counts = {}
        for g in genotypes.values():
            status_counts[g['status']] = status_counts.get(g['status'], 0) + 1
        print(f"Status summary: {status_counts}")

        # Generate output
        print(f"Writing output to {args.output}...")
        write_output(genotypes, args.output, args.format)

        # Print fingerprint
        fingerprint = generate_fingerprint(genotypes)
        fingerprint_hash = generate_hash_fingerprint(genotypes)
        print(f"\nSample Fingerprint: {fingerprint}")
        print(f"Fingerprint Hash: {fingerprint_hash}")

        # Cleanup
        if not args.keep_vcf:
            os.remove(vcf_file)
            print(f"Removed temporary VCF file")

    finally:
        os.remove(positions_file)


if __name__ == "__main__":
    main()