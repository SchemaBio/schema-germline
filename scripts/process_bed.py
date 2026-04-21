#!/usr/bin/env python3
"""
Process BED file: normalize chromosome names and fix zero-length regions.

Functions:
1. Remove 'chr' prefix from chromosome names (e.g., chr1 -> 1)
2. Convert 'chrM' to 'MT'
3. Fix regions where end == start by incrementing end by 1
"""

import argparse
import sys
from pathlib import Path


def process_chromosome(chrom: str) -> str:
    """Normalize chromosome name by removing 'chr' prefix."""
    if chrom == "chrM":
        return "MT"
    elif chrom.startswith("chr"):
        return chrom[3:]
    return chrom


def process_bed_file(input_file: str, output_file: str) -> dict:
    """
    Process BED file and write normalized output.

    Args:
        input_file: Path to input BED file
        output_file: Path to output BED file

    Returns:
        Dictionary with processing statistics
    """
    stats = {
        "total_lines": 0,
        "chr_removed": 0,
        "chrM_to_MT": 0,
        "zero_length_fixed": 0
    }

    with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
        for line in f_in:
            line = line.rstrip("\n")
            stats["total_lines"] += 1

            # Skip empty lines and comments
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                f_out.write(line + "\n")
                continue

            fields = line.split("\t")
            if len(fields) < 3:
                f_out.write(line + "\n")
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            # Process chromosome name
            original_chrom = chrom
            new_chrom = process_chromosome(chrom)

            if original_chrom != new_chrom:
                fields[0] = new_chrom
                if original_chrom == "chrM":
                    stats["chrM_to_MT"] += 1
                else:
                    stats["chr_removed"] += 1

            # Fix zero-length regions (end == start)
            if end == start:
                fields[2] = str(end + 1)
                stats["zero_length_fixed"] += 1

            f_out.write("\t".join(fields) + "\n")

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Process BED file: normalize chromosome names and fix zero-length regions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.bed output.bed
  %(prog)s -i targets.bed -o targets_normalized.bed
        """
    )
    parser.add_argument("-i", "--input", required=True, help="Input BED file")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")

    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)

    # Create output directory if needed
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Processing BED file: {args.input}")

    stats = process_bed_file(args.input, args.output)

    print(f"\nProcessing complete:")
    print(f"  Total lines processed: {stats['total_lines']}")
    print(f"  'chr' prefix removed: {stats['chr_removed']}")
    print(f"  'chrM' -> 'MT' converted: {stats['chrM_to_MT']}")
    print(f"  Zero-length regions fixed: {stats['zero_length_fixed']}")
    print(f"\nOutput written to: {args.output}")


if __name__ == "__main__":
    main()