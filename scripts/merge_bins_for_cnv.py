#!/usr/bin/env python3
"""
Bin merging script for CNVkit .cnr files.

Merges small exon-level bins into fixed-size bins (e.g., 10kb, 20kb, 50kb)
before segmentation for detecting large CNV events.

Usage:
    python merge_bins_for_cnv.py -i input.cnr -o merged.cnr --bin-size 20000

The merged file can then be used with gene_level_segment.py for segmentation.
"""

import argparse
import sys
import statistics
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass


@dataclass
class Bin:
    """Original bin from .cnr file."""
    chromosome: str
    start: int
    end: int
    gene: str
    log2: float
    depth: float
    weight: float


@dataclass
class MergedBin:
    """Merged bin from multiple original bins."""
    chromosome: str
    start: int
    end: int
    genes: list  # List of gene names
    log2: float  # Weighted average
    depth: float  # Weighted average by weight
    weight: float  # Sum of weights
    bin_count: int  # Number of original bins merged
    cn: float  # Copy number, 2 decimal places
    cnv_call: str  # DUP, DEL, or Normal based on thresholds


def calculate_cn(log2: float) -> float:
    """Calculate copy number from log2 ratio. CN = 2^(log2 + 1) = 2 * 2^log2."""
    import math
    cn = 2 ** (log2 + 1)
    return round(cn, 2)


def call_cnv(cn: float, dup_threshold: float, del_threshold: float) -> str:
    """
    Determine CNV call based on copy number and thresholds.

    Args:
        cn: Copy number value
        dup_threshold: CN threshold for DUP (CN >= this value)
        del_threshold: CN threshold for DEL (CN <= this value)

    Returns:
        "DUP", "DEL", or "Normal"
    """
    if cn >= dup_threshold:
        return "DUP"
    elif cn <= del_threshold:
        return "DEL"
    else:
        return "Normal"


def read_cnr_file(input_file: str, skip_antitarget: bool = True) -> list[Bin]:
    """Read CNVkit .cnr file."""
    bins = []
    with open(input_file, "r") as f:
        header = f.readline()
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 7:
                continue

            gene_field = fields[3]
            # Optionally skip Antitarget regions
            if skip_antitarget and gene_field == "Antitarget":
                continue

            bin_data = Bin(
                chromosome=fields[0],
                start=int(fields[1]),
                end=int(fields[2]),
                gene=gene_field,
                log2=float(fields[4]),
                depth=float(fields[5]),
                weight=float(fields[6])
            )
            bins.append(bin_data)

    return bins


def merge_bins_by_size(bins: list[Bin], bin_size: int, dup_threshold: float = 2.5, del_threshold: float = 1.5) -> list[MergedBin]:
    """
    Merge bins into fixed-size bins by genomic position.

    Args:
        bins: List of original Bin objects
        bin_size: Target bin size in base pairs (e.g., 10000, 20000, 50000)
        dup_threshold: CN threshold for DUP call (default: 2.5)
        del_threshold: CN threshold for DEL call (default: 1.5)

    Returns:
        List of MergedBin objects
    """
    # Group bins by chromosome
    chrom_bins = defaultdict(list)
    for bin_data in bins:
        chrom_bins[bin_data.chromosome].append(bin_data)

    # Sort bins within each chromosome by start position
    for chrom in chrom_bins:
        chrom_bins[chrom].sort(key=lambda b: b.start)

    merged_bins = []

    for chrom, chrom_bin_list in chrom_bins.items():
        if not chrom_bin_list:
            continue

        # Assign each bin to a window based on its center position
        # This ensures each bin is only assigned to one window
        window_bins = defaultdict(list)

        for b in chrom_bin_list:
            # Use bin center to determine which window it belongs to
            bin_center = (b.start + b.end) // 2
            window_idx = bin_center // bin_size
            window_bins[window_idx].append(b)

        # Merge bins within each window
        for window_idx, window_bin_list in sorted(window_bins.items()):
            if not window_bin_list:
                continue

            # Calculate weighted statistics
            total_weight = sum(b.weight for b in window_bin_list)
            if total_weight > 0:
                log2_weighted = sum(b.log2 * b.weight for b in window_bin_list) / total_weight
                depth_weighted = sum(b.depth * b.weight for b in window_bin_list) / total_weight
            else:
                log2_weighted = statistics.mean(b.log2 for b in window_bin_list)
                depth_weighted = statistics.mean(b.depth for b in window_bin_list)

            # Collect gene names (unique)
            genes = list(set(b.gene for b in window_bin_list if b.gene != "Antitarget"))

            # Window boundaries
            win_start = window_idx * bin_size
            win_end = (window_idx + 1) * bin_size

            # Actual start/end based on bins (clipped to window for cleaner boundaries)
            actual_start = max(min(b.start for b in window_bin_list), win_start)
            actual_end = min(max(b.end for b in window_bin_list), win_end)

            # Calculate CN and CNV call
            cn = calculate_cn(log2_weighted)
            cnv_call = call_cnv(cn, dup_threshold, del_threshold)

            merged = MergedBin(
                chromosome=chrom,
                start=actual_start,
                end=actual_end,
                genes=genes,
                log2=log2_weighted,
                depth=depth_weighted,
                weight=total_weight,
                bin_count=len(window_bin_list),
                cn=cn,
                cnv_call=cnv_call
            )
            merged_bins.append(merged)

    # Sort by chromosome and position
    def sort_key(mb: MergedBin) -> tuple:
        chrom = mb.chromosome
        if chrom.isdigit():
            return (int(chrom), mb.start)
        elif chrom == "X":
            return (23, mb.start)
        elif chrom == "Y":
            return (24, mb.start)
        elif chrom == "MT":
            return (25, mb.start)
        else:
            return (26, mb.start)

    merged_bins.sort(key=sort_key)

    return merged_bins


def merge_bins_by_count(bins: list[Bin], bin_count: int, dup_threshold: float = 2.5, del_threshold: float = 1.5) -> list[MergedBin]:
    """
    Merge adjacent bins by fixed count (e.g., every 5 bins merged into 1).

    Args:
        bins: List of original Bin objects
        bin_count: Number of bins to merge together
        dup_threshold: CN threshold for DUP call (default: 2.5)
        del_threshold: CN threshold for DEL call (default: 1.5)

    Returns:
        List of MergedBin objects
    """
    # Group bins by chromosome
    chrom_bins = defaultdict(list)
    for bin_data in bins:
        chrom_bins[bin_data.chromosome].append(bin_data)

    # Sort bins within each chromosome by start position
    for chrom in chrom_bins:
        chrom_bins[chrom].sort(key=lambda b: b.start)

    merged_bins = []

    for chrom, chrom_bin_list in chrom_bins.items():
        if not chrom_bin_list:
            continue

        # Merge every N bins
        for i in range(0, len(chrom_bin_list), bin_count):
            chunk = chrom_bin_list[i:i + bin_count]
            if not chunk:
                continue

            # Calculate weighted statistics
            total_weight = sum(b.weight for b in chunk)
            if total_weight > 0:
                log2_weighted = sum(b.log2 * b.weight for b in chunk) / total_weight
                depth_weighted = sum(b.depth * b.weight for b in chunk) / total_weight
            else:
                log2_weighted = statistics.mean(b.log2 for b in chunk)
                depth_weighted = statistics.mean(b.depth for b in chunk)

            genes = list(set(b.gene for b in chunk if b.gene != "Antitarget"))

            # Calculate CN and CNV call
            cn = calculate_cn(log2_weighted)
            cnv_call = call_cnv(cn, dup_threshold, del_threshold)

            merged = MergedBin(
                chromosome=chrom,
                start=min(b.start for b in chunk),
                end=max(b.end for b in chunk),
                genes=genes,
                log2=log2_weighted,
                depth=depth_weighted,
                weight=total_weight,
                bin_count=len(chunk),
                cn=cn,
                cnv_call=cnv_call
            )
            merged_bins.append(merged)

    # Sort
    def sort_key(mb: MergedBin) -> tuple:
        chrom = mb.chromosome
        if chrom.isdigit():
            return (int(chrom), mb.start)
        elif chrom == "X":
            return (23, mb.start)
        elif chrom == "Y":
            return (24, mb.start)
        elif chrom == "MT":
            return (25, mb.start)
        else:
            return (26, mb.start)

    merged_bins.sort(key=sort_key)

    return merged_bins


def write_merged_cnr(merged_bins: list[MergedBin], output_file: str):
    """Write merged bins to CNVkit .cnr format with CN and CNV_call columns."""
    with open(output_file, "w") as f:
        # Header (CNV_call as 4th column)
        header = ["chromosome", "start", "end", "CNV_call", "log2", "depth", "weight", "CN"]
        f.write("\t".join(header) + "\n")

        for mb in merged_bins:
            row = [
                mb.chromosome,
                str(mb.start),
                str(mb.end),
                mb.cnv_call,
                f"{mb.log2:.6f}",
                f"{mb.depth:.2f}",
                f"{mb.weight:.6f}",
                f"{mb.cn:.2f}"
            ]
            f.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Merge CNVkit .cnr bins into fixed-size bins for large CNV detection.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Merge into 20kb bins (recommended for large CNV detection)
  %(prog)s -i proband.cnvkit.cnr -o proband_20kb.cnr --bin-size 20000

  # Merge into 50kb bins for very large CNVs
  %(prog)s -i proband.cnvkit.cnr -o proband_50kb.cnr --bin-size 50000

  # Merge every 5 bins together
  %(prog)s -i proband.cnvkit.cnr -o proband_merged.cnr --bin-count 5

  # Include antitarget regions
  %(prog)s -i proband.cnvkit.cnr -o merged.cnr --bin-size 20000 --keep-antitarget

  # Custom CNV thresholds (CN >= 3.0 for DUP, CN <= 1.0 for DEL)
  %(prog)s -i proband.cnvkit.cnr -o merged.cnr --bin-size 20000 --dup-threshold 3.0 --del-threshold 1.0

Recommended bin sizes:
  - 10kb:  Good for detecting medium CNVs (~50kb-1Mb)
  - 20kb:  Good for large CNVs (~100kb-10Mb) [recommended]
  - 50kb:  Good for very large CNVs (>1Mb), reduces noise but loses resolution

CNV calling:
  - CN (Copy Number) is calculated as: CN = 2^(log2 + 1) = 2 * 2^log2
  - DUP: CN >= dup_threshold (default: 2.5)
  - DEL: CN <= del_threshold (default: 1.5)
  - Normal: del_threshold < CN < dup_threshold
        """
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Input CNVkit .cnr file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output merged .cnr file")

    # Merging method
    method_group = parser.add_mutually_exclusive_group(required=True)
    method_group.add_argument("--bin-size", type=int,
                              help="Fixed bin size in base pairs (e.g., 10000, 20000, 50000)")
    method_group.add_argument("--bin-count", type=int,
                              help="Number of adjacent bins to merge (e.g., 5)")

    parser.add_argument("--keep-antitarget", action="store_true",
                        help="Keep Antitarget regions (default: skip)")

    # CNV thresholds
    parser.add_argument("--dup-threshold", type=float, default=2.5,
                        help="CN threshold for DUP call (default: 2.5, CN >= this value)")
    parser.add_argument("--del-threshold", type=float, default=1.5,
                        help="CN threshold for DEL call (default: 1.5, CN <= this value)")

    args = parser.parse_args()

    # Read input
    print(f"Reading CNVkit .cnr file: {args.input}")
    bins = read_cnr_file(args.input, skip_antitarget=not args.keep_antitarget)
    print(f"  Total bins: {len(bins)}")

    # Merge bins
    print(f"\nMerging bins...")
    print(f"  DUP threshold: CN >= {args.dup_threshold}")
    print(f"  DEL threshold: CN <= {args.del_threshold}")
    if args.bin_size:
        merged_bins = merge_bins_by_size(bins, args.bin_size, args.dup_threshold, args.del_threshold)
        print(f"  Method: Fixed size ({args.bin_size} bp)")
    else:
        merged_bins = merge_bins_by_count(bins, args.bin_count, args.dup_threshold, args.del_threshold)
        print(f"  Method: Fixed count ({args.bin_count} bins)")

    print(f"  Merged bins: {len(merged_bins)}")

    # Statistics
    avg_original_size = statistics.mean(b.end - b.start for b in bins)
    avg_merged_size = statistics.mean(mb.end - mb.start for mb in merged_bins)
    avg_bins_merged = statistics.mean(mb.bin_count for mb in merged_bins)

    print(f"\nStatistics:")
    print(f"  Original avg bin size: {avg_original_size:.0f} bp")
    print(f"  Merged avg bin size: {avg_merged_size:.0f} bp")
    print(f"  Avg bins merged per window: {avg_bins_merged:.1f}")

    # CNV statistics
    dup_count = sum(1 for mb in merged_bins if mb.cnv_call == "DUP")
    del_count = sum(1 for mb in merged_bins if mb.cnv_call == "DEL")
    normal_count = sum(1 for mb in merged_bins if mb.cnv_call == "Normal")
    print(f"\nCNV Calls:")
    print(f"  DUP: {dup_count} ({100*dup_count/len(merged_bins):.1f}%)")
    print(f"  DEL: {del_count} ({100*del_count/len(merged_bins):.1f}%)")
    print(f"  Normal: {normal_count} ({100*normal_count/len(merged_bins):.1f}%)")

    # Write output
    print(f"\nWriting merged file: {args.output}")
    write_merged_cnr(merged_bins, args.output)

    print(f"\nDone! Use the merged file with gene_level_segment.py for CNV detection.")


if __name__ == "__main__":
    main()