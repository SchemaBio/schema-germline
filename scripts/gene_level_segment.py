#!/usr/bin/env python3
"""
Gene-level CNV segmentation for CNVkit .cnr files.

This script implements CBS (Circular Binary Segmentation) algorithm followed by
gene-level aggregation for germline whole-exome CNV detection.

Algorithm:
1. Optional bin merging: Merge small bins into fixed-size bins (e.g., 10k, 20k)
2. CBS segmentation: Identify copy number change points on each chromosome
3. Gene mapping: Map segments to genes based on coverage overlap
4. CNV classification: Classify as Normal/Deletion/Duplication with confidence

Output includes:
- Gene-level CNV calls with segment-level statistics
- CNV type classification (DEL/DUP/NEUTRAL)
- Confidence metrics (p-value, segment quality)
"""

import argparse
import sys
import math
import statistics
from pathlib import Path
from collections import defaultdict, Counter
from typing import NamedTuple, Optional, Union
from dataclasses import dataclass


# ============================================================================
# Data structures
# ============================================================================

@dataclass
class Bin:
    """Bin-level data from .cnr file."""
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
    genes: list  # List of original gene annotations
    log2: float  # Weighted average
    depth: float  # Weighted average
    weight: float  # Sum of weights
    bin_count: int  # Number of original bins merged
    original_bins: list  # List of original Bin objects


@dataclass
class Segment:
    """CBS segment result."""
    chromosome: str
    start: int
    end: int
    log2_mean: float
    log2_std: float
    bin_count: int
    weight_sum: float
    probes: list  # List of bin indices


class GeneCNV(NamedTuple):
    """Gene-level CNV result."""
    chromosome: str
    gene: str
    transcript: str
    ensembl: str
    exon_count: int
    start: int
    end: int
    log2_mean: float
    log2_median: float
    log2_std: float
    cn: float  # Estimated copy number (2 decimal places)
    depth_mean: float
    weight_mean: float
    bin_count: int
    segment_count: int
    coverage_ratio: float  # Ratio of gene covered by segments
    cnv_call: str  # DEL, DUP, Normal (based on thresholds)
    confidence: str  # HIGH, MEDIUM, LOW
    p_value: float


# ============================================================================
# CBS Algorithm Implementation
# ============================================================================

def calculate_cusum(values: list[float], mean: float) -> list[float]:
    """Calculate cumulative sum of deviations from mean."""
    cusum = []
    cum = 0.0
    for v in values:
        cum += (v - mean)
        cusum.append(cum)
    return cusum


def find_max_cusum_point(cusum: list[float]) -> tuple[int, float]:
    """Find point with maximum absolute cusum deviation."""
    max_idx = 0
    max_val = 0.0
    for i, c in enumerate(cusum):
        if abs(c) > abs(max_val):
            max_val = c
            max_idx = i
    return max_idx, max_val


def calculate_cbs_statistic(values: list[float], k: int) -> float:
    """
    Calculate CBS test statistic for a potential change point at position k.

    The CBS statistic tests if there's a significant difference in mean
    between the regions before and after position k.
    """
    n = len(values)
    if k <= 0 or k >= n - 1:
        return 0.0

    # Calculate means before and after k
    mean1 = statistics.mean(values[:k+1])
    mean2 = statistics.mean(values[k+1:])

    # Calculate variance (pooled)
    if n > 2:
        var1 = statistics.variance(values[:k+1]) if k > 0 else 0.0
        var2 = statistics.variance(values[k+1:]) if k < n-1 else 0.0

        # Pooled variance estimate
        n1 = k + 1
        n2 = n - k - 1
        pooled_var = (n1 * var1 + n2 * var2) / (n1 + n2 - 2) if n1 + n2 > 2 else 0.0

        # t-statistic
        if pooled_var > 0:
            se = math.sqrt(pooled_var * (1/n1 + 1/n2))
            t_stat = abs(mean1 - mean2) / se
            return t_stat
    return abs(mean1 - mean2)


def calculate_p_value(t_stat: float, n: int, k: int) -> float:
    """
    Approximate p-value for CBS test statistic.

    Uses t-distribution approximation with Bonferroni correction for
    multiple testing (testing all possible change points).
    """
    if n <= 2 or t_stat <= 0:
        return 1.0

    # Degrees of freedom
    df = n - 2

    # Simple approximation of p-value using t-distribution
    # For large n, t-distribution approaches normal
    if df > 30:
        # Normal approximation
        p = 2 * (1 - normal_cdf(abs(t_stat)))
    else:
        # t-distribution approximation (simplified)
        p = t_distribution_p_value_approx(t_stat, df)

    # Bonferroni correction for testing n-1 possible change points
    p_corrected = min(p * (n - 1), 1.0)

    return p_corrected


def normal_cdf(x: float) -> float:
    """Approximate standard normal CDF."""
    # Abramowitz and Stegun approximation
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    sign = 1 if x >= 0 else -1
    x = abs(x) / math.sqrt(2)

    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-x * x)

    return 0.5 * (1.0 + sign * y)


def t_distribution_p_value_approx(t: float, df: int) -> float:
    """
    Approximate two-tailed p-value for t-distribution.

    Uses a simple approximation based on the relationship between
    t-distribution and normal distribution.
    """
    t = abs(t)

    # For large df, use normal approximation
    if df > 60:
        return 2 * (1 - normal_cdf(t))

    # For moderate df, use adjusted normal approximation
    # This adjusts for the heavier tails of t-distribution
    adjusted_t = t * math.sqrt(df / (df - 2)) if df > 2 else t
    return 2 * (1 - normal_cdf(adjusted_t))


def cbs_segment_recursive(
    values: list[float],
    positions: list[int],
    weights: list[float],
    alpha: float = 0.01,
    min_bins: int = 3,
    max_depth: int = 10,
    depth: int = 0
) -> list[tuple[int, int]]:
    """
    Recursive CBS segmentation.

    Args:
        values: Log2 values for bins
        positions: Genomic positions (midpoints)
        weights: Weights for each bin
        alpha: Significance threshold for change point
        min_bins: Minimum bins per segment
        max_depth: Maximum recursion depth
        depth: Current recursion depth

    Returns:
        List of (start_idx, end_idx) segment boundaries
    """
    n = len(values)

    # Stop conditions
    if n < min_bins * 2 or depth >= max_depth:
        return [(0, n - 1)]

    # Calculate overall mean
    total_weight = sum(weights)
    if total_weight <= 0:
        weights = [1.0] * n
        total_weight = n

    mean = sum(v * w for v, w in zip(values, weights)) / total_weight

    # Calculate cusum
    cusum = calculate_cusum(values, mean)

    # Find maximum cusum point
    k, max_cusum = find_max_cusum_point(cusum)

    # Calculate test statistic
    t_stat = calculate_cbs_statistic(values, k)

    # Calculate p-value
    p_value = calculate_p_value(t_stat, n, k)

    # Test if change point is significant
    if p_value > alpha or k < min_bins or k > n - min_bins - 1:
        return [(0, n - 1)]

    # Recursively segment left and right regions
    left_segments = cbs_segment_recursive(
        values[:k+1],
        positions[:k+1],
        weights[:k+1],
        alpha, min_bins, max_depth, depth + 1
    )

    right_segments = cbs_segment_recursive(
        values[k+1:],
        positions[k+1:],
        weights[k+1:],
        alpha, min_bins, max_depth, depth + 1
    )

    # Combine segments
    segments = []
    for start, end in left_segments:
        segments.append((start, end))
    for start, end in right_segments:
        segments.append((start + k + 1, end + k + 1))

    return segments


def run_cbs_on_chromosome(
    bins: list[Bin],
    alpha: float = 0.01,
    min_bins: int = 3
) -> list[Segment]:
    """
    Run CBS segmentation on bins from one chromosome.

    Args:
        bins: List of Bin objects for one chromosome
        alpha: Significance threshold
        min_bins: Minimum bins per segment

    Returns:
        List of Segment objects
    """
    if len(bins) < min_bins:
        # Return single segment if too few bins
        log2_values = [b.log2 for b in bins]
        weights = [b.weight for b in bins]
        total_weight = sum(weights)
        if total_weight > 0:
            mean = sum(l * w for l, w in zip(log2_values, weights)) / total_weight
        else:
            mean = statistics.mean(log2_values) if log2_values else 0.0

        std = statistics.stdev(log2_values) if len(log2_values) > 1 else 0.0

        return [Segment(
            chromosome=bins[0].chromosome if bins else "",
            start=min(b.start for b in bins) if bins else 0,
            end=max(b.end for b in bins) if bins else 0,
            log2_mean=mean,
            log2_std=std,
            bin_count=len(bins),
            weight_sum=sum(b.weight for b in bins),
            probes=list(range(len(bins)))
        )]

    # Prepare data
    log2_values = [b.log2 for b in bins]
    positions = [(b.start + b.end) // 2 for b in bins]  # Midpoints
    weights = [b.weight for b in bins]

    # Run CBS
    segment_boundaries = cbs_segment_recursive(
        log2_values, positions, weights, alpha, min_bins
    )

    # Build segment objects
    segments = []
    for start_idx, end_idx in segment_boundaries:
        segment_bins = bins[start_idx:end_idx + 1]
        seg_log2 = [b.log2 for b in segment_bins]
        seg_weights = [b.weight for b in segment_bins]

        total_weight = sum(seg_weights)
        if total_weight > 0:
            mean = sum(l * w for l, w in zip(seg_log2, seg_weights)) / total_weight
        else:
            mean = statistics.mean(seg_log2)

        std = statistics.stdev(seg_log2) if len(seg_log2) > 1 else 0.0

        segment = Segment(
            chromosome=segment_bins[0].chromosome,
            start=min(b.start for b in segment_bins),
            end=max(b.end for b in segment_bins),
            log2_mean=mean,
            log2_std=std,
            bin_count=len(segment_bins),
            weight_sum=total_weight,
            probes=list(range(start_idx, end_idx + 1))
        )
        segments.append(segment)

    return segments


# ============================================================================
# Gene parsing and mapping
# ============================================================================

def parse_gene_field(gene_field: str) -> tuple[str, str, str, str]:
    """
    Parse the gene annotation field from CNVkit .cnr file.

    Format: GeneName|TranscriptID|ENSEMBL_ID|ExonNumber|Strand|Cytoband

    Returns:
        Tuple of (gene_name, transcript, ensembl, cytoband)
    """
    parts = gene_field.split("|")
    if len(parts) >= 6:
        return parts[0], parts[1], parts[2], parts[5]
    elif len(parts) >= 4:
        return parts[0], parts[1], parts[2], ""
    elif len(parts) >= 1:
        return parts[0], "", "", ""
    return gene_field, "", "", ""


def read_cnr_file(input_file: str) -> list[Bin]:
    """Read CNVkit .cnr file."""
    bins = []
    with open(input_file, "r") as f:
        header = f.readline()
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 7:
                continue

            gene_field = fields[3]
            # Skip Antitarget regions
            if gene_field == "Antitarget":
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


def calculate_overlap_ratio(seg_start: int, seg_end: int,
                           gene_start: int, gene_end: int) -> float:
    """Calculate overlap ratio between segment and gene."""
    overlap_start = max(seg_start, gene_start)
    overlap_end = min(seg_end, gene_end)

    if overlap_start >= overlap_end:
        return 0.0

    overlap_length = overlap_end - overlap_start
    gene_length = gene_end - gene_start

    return overlap_length / gene_length if gene_length > 0 else 0.0


def calculate_cn(log2_mean: float, ploidy: int = 2) -> float:
    """
    Calculate copy number from log2 value.

    Formula: CN = ploidy * 2^log2

    For diploid (ploidy=2):
    - CN=0: log2 = log(0/2) = -∞ (practically very negative)
    - CN=1: log2 = log(1/2) ≈ -0.5
    - CN=2: log2 = log(2/2) = 0
    - CN=3: log2 = log(3/2) ≈ 0.58
    - CN=4: log2 = log(4/2) ≈ 1.0

    Returns:
        Copy number as float (rounded to 2 decimal places)
    """
    cn = ploidy * math.pow(2, log2_mean)
    return max(0.0, round(cn, 2))


def classify_cnv(
    log2_mean: float,
    bin_count: int,
    log2_std: float,
    del_threshold: float,
    dup_threshold: float
) -> tuple[str, str, float]:
    """
    Classify CNV type based on log2 value and thresholds.

    Args:
        log2_mean: Mean log2 ratio
        bin_count: Number of bins covering the gene
        log2_std: Standard deviation of log2 values
        del_threshold: Log2 threshold for calling deletion (negative)
        dup_threshold: Log2 threshold for calling duplication (positive)

    Returns:
        Tuple of (cnv_call, confidence, cn)

    cnv_call:
    - DEL: log2 < del_threshold
    - DUP: log2 > dup_threshold
    - Normal: otherwise

    Confidence levels:
    - HIGH: log2 strongly deviates, low std, many bins
    - MEDIUM: moderate deviation
    - LOW: borderline values, high variance, few bins
    """
    # Calculate copy number
    cn = calculate_cn(log2_mean)

    # Base classification based on thresholds
    if log2_mean < del_threshold:
        cnv_call = "DEL"
    elif log2_mean > dup_threshold:
        cnv_call = "DUP"
    else:
        cnv_call = "Normal"

    # Confidence assessment
    confidence_score = 0

    # Distance from threshold contribution
    if cnv_call == "DEL":
        distance = abs(log2_mean - del_threshold)
        if distance > 0.3:
            confidence_score += 2
        elif distance > 0.15:
            confidence_score += 1
    elif cnv_call == "DUP":
        distance = abs(log2_mean - dup_threshold)
        if distance > 0.2:
            confidence_score += 2
        elif distance > 0.1:
            confidence_score += 1
    else:
        # Normal: high confidence if clearly in normal range
        min_dist = min(abs(log2_mean - del_threshold), abs(log2_mean - dup_threshold))
        if min_dist > 0.2:
            confidence_score += 2
        elif min_dist > 0.1:
            confidence_score += 1

    # Variance contribution (low variance = higher confidence)
    if log2_std < 0.15:
        confidence_score += 1
    elif log2_std > 0.4:
        confidence_score -= 1

    # Bin count contribution (more bins = higher confidence)
    if bin_count >= 5:
        confidence_score += 1
    elif bin_count < 3:
        confidence_score -= 1

    # Final confidence level
    if confidence_score >= 3:
        confidence = "HIGH"
    elif confidence_score >= 1:
        confidence = "MEDIUM"
    else:
        confidence = "LOW"

    return cnv_call, confidence, cn


def segment_to_gene_cnv(
    segments: list[Segment],
    bins: list[Bin],
    chromosome: str,
    del_threshold: float = -0.4,
    dup_threshold: float = 0.3
) -> list[GeneCNV]:
    """
    Map CBS segments to genes and calculate gene-level CNV statistics.

    For each gene, find overlapping segments and aggregate:
    - If one segment covers most of gene -> use segment directly
    - If multiple segments overlap -> use weighted average based on coverage

    Args:
        segments: List of CBS segments
        bins: List of Bin objects
        chromosome: Chromosome being processed
        del_threshold: Log2 threshold for calling deletion
        dup_threshold: Log2 threshold for calling duplication
    """
    # Group bins by gene
    gene_bins = defaultdict(list)
    for i, bin_data in enumerate(bins):
        gene_name, transcript, ensembl, cytoband = parse_gene_field(bin_data.gene)
        gene_bins[gene_name].append({
            "bin": bin_data,
            "index": i,
            "transcript": transcript,
            "ensembl": ensembl,
            "cytoband": cytoband
        })

    results = []

    for gene_name, gene_data_list in gene_bins.items():
        gene_bins_list = [d["bin"] for d in gene_data_list]
        gene_start = min(b.start for b in gene_bins_list)
        gene_end = max(b.end for b in gene_bins_list)

        # Find overlapping segments
        overlapping_segments = []
        for seg in segments:
            if seg.chromosome == chromosome:
                overlap_ratio = calculate_overlap_ratio(seg.start, seg.end, gene_start, gene_end)
                if overlap_ratio > 0.1:  # At least 10% overlap
                    overlapping_segments.append((seg, overlap_ratio))

        if not overlapping_segments:
            continue

        # Calculate gene statistics
        # Use the segment with maximum overlap as primary
        # But also aggregate all overlapping bins for robust statistics

        # Get primary segment
        primary_segment = max(overlapping_segments, key=lambda x: x[1])
        seg, coverage_ratio = primary_segment

        # Aggregate all gene bins for robust statistics
        log2_values = [b.log2 for b in gene_bins_list]
        weights = [b.weight for b in gene_bins_list]

        total_weight = sum(weights)
        if total_weight > 0:
            log2_weighted_mean = sum(l * w for l, w in zip(log2_values, weights)) / total_weight
        else:
            log2_weighted_mean = statistics.mean(log2_values)

        log2_median = statistics.median(log2_values)
        log2_std = statistics.stdev(log2_values) if len(log2_values) > 1 else 0.0
        depth_mean = statistics.mean(b.depth for b in gene_bins_list)
        weight_mean = statistics.mean(weights)

        # Get transcript/ensembl (most common)
        transcripts = [d["transcript"] for d in gene_data_list if d["transcript"]]
        ensembls = [d["ensembl"] for d in gene_data_list if d["ensembl"]]

        transcript = Counter(transcripts).most_common(1)[0][0] if transcripts else ""
        ensembl = Counter(ensembls).most_common(1)[0][0] if ensembls else ""

        # Count exons
        exon_count = len(set(
            b.gene.split("|")[3] for b in gene_bins_list
            if len(b.gene.split("|")) >= 4
        ))
        if exon_count == 0:
            exon_count = len(gene_bins_list)

        # Calculate p-value for CNV call
        n = len(log2_values)
        if n > 1 and log2_std > 0:
            t_stat = abs(log2_weighted_mean) / (log2_std / math.sqrt(n))
            df = n - 1
            if df > 30:
                p_value = 2 * (1 - normal_cdf(t_stat))
            else:
                p_value = t_distribution_p_value_approx(t_stat, df)
            p_value = min(p_value * n, 1.0)
        else:
            p_value = 1.0

        # Classify CNV with thresholds
        cnv_call, confidence, cn = classify_cnv(
            log2_weighted_mean, len(gene_bins_list), log2_std,
            del_threshold, dup_threshold
        )

        gene_cnv = GeneCNV(
            chromosome=chromosome,
            gene=gene_name,
            transcript=transcript,
            ensembl=ensembl,
            exon_count=exon_count,
            start=gene_start,
            end=gene_end,
            log2_mean=log2_weighted_mean,
            log2_median=log2_median,
            log2_std=log2_std,
            cn=cn,
            depth_mean=depth_mean,
            weight_mean=weight_mean,
            bin_count=len(gene_bins_list),
            segment_count=len(overlapping_segments),
            coverage_ratio=coverage_ratio,
            cnv_call=cnv_call,
            confidence=confidence,
            p_value=p_value
        )

        results.append(gene_cnv)

    return results


# ============================================================================
# Main workflow
# ============================================================================

def run_cbs_segmentation(
    bins: list[Bin],
    alpha: float = 0.01,
    min_bins: int = 3
) -> dict[str, list[Segment]]:
    """
    Run CBS segmentation on all chromosomes.

    Returns:
        Dictionary mapping chromosome to list of segments
    """
    # Group bins by chromosome
    chrom_bins = defaultdict(list)
    for bin_data in bins:
        chrom_bins[bin_data.chromosome].append(bin_data)

    # Run CBS on each chromosome
    all_segments = {}
    for chrom, chrom_bin_list in chrom_bins.items():
        print(f"  Processing chromosome {chrom}: {len(chrom_bin_list)} bins")
        segments = run_cbs_on_chromosome(chrom_bin_list, alpha, min_bins)
        all_segments[chrom] = segments
        print(f"    -> {len(segments)} segments")

    return all_segments


def gene_level_cnv(
    bins: list[Bin],
    segments: dict[str, list[Segment]],
    del_threshold: float = -0.4,
    dup_threshold: float = 0.3
) -> list[GeneCNV]:
    """
    Map segments to genes and generate gene-level CNV calls.

    Args:
        bins: List of Bin objects
        segments: Dictionary mapping chromosome to segments
        del_threshold: Log2 threshold for calling deletion
        dup_threshold: Log2 threshold for calling duplication
    """
    # Group bins by chromosome
    chrom_bins = defaultdict(list)
    for bin_data in bins:
        chrom_bins[bin_data.chromosome].append(bin_data)

    results = []
    for chrom, chrom_bin_list in chrom_bins.items():
        chrom_segments = segments.get(chrom, [])
        gene_cnvs = segment_to_gene_cnv(
            chrom_segments, chrom_bin_list, chrom,
            del_threshold, dup_threshold
        )
        results.extend(gene_cnvs)

    # Sort by chromosome and position
    def sort_key(cnv: GeneCNV) -> tuple:
        chrom = cnv.chromosome
        if chrom.isdigit():
            return (int(chrom), cnv.start)
        elif chrom == "X":
            return (23, cnv.start)
        elif chrom == "Y":
            return (24, cnv.start)
        elif chrom == "MT":
            return (25, cnv.start)
        else:
            return (26, cnv.start)

    results.sort(key=sort_key)

    return results


def write_output(cnvs: list[GeneCNV], output_file: str):
    """Write gene-level CNV results to file."""
    with open(output_file, "w") as f:
        header = [
            "chromosome", "start", "end", "cnv_call", "gene", "transcript", "ensembl", "exon_count",
            "log2_mean", "log2_median", "log2_std", "cn",
            "depth_mean", "weight_mean", "bin_count", "segment_count",
            "coverage_ratio", "confidence", "p_value"
        ]
        f.write("\t".join(header) + "\n")

        for cnv in cnvs:
            row = [
                cnv.chromosome,
                str(cnv.start),
                str(cnv.end),
                cnv.cnv_call,
                cnv.gene,
                cnv.transcript,
                cnv.ensembl,
                str(cnv.exon_count),
                f"{cnv.log2_mean:.6f}",
                f"{cnv.log2_median:.6f}",
                f"{cnv.log2_std:.6f}",
                f"{cnv.cn:.2f}",
                f"{cnv.depth_mean:.2f}",
                f"{cnv.weight_mean:.6f}",
                str(cnv.bin_count),
                str(cnv.segment_count),
                f"{cnv.coverage_ratio:.4f}",
                cnv.confidence,
                f"{cnv.p_value:.6f}"
            ]
            f.write("\t".join(row) + "\n")


def write_segments_file(segments: dict[str, list[Segment]], output_file: str):
    """Write CBS segments to file."""
    with open(output_file, "w") as f:
        header = ["chromosome", "start", "end", "log2_mean", "log2_std", "bin_count", "weight_sum"]
        f.write("\t".join(header) + "\n")

        for chrom in sorted(segments.keys(), key=lambda x: int(x) if x.isdigit() else 99):
            for seg in segments[chrom]:
                row = [
                    seg.chromosome,
                    str(seg.start),
                    str(seg.end),
                    f"{seg.log2_mean:.6f}",
                    f"{seg.log2_std:.6f}",
                    str(seg.bin_count),
                    f"{seg.weight_sum:.6f}"
                ]
                f.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="CBS-based gene-level CNV segmentation for CNVkit .cnr files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i proband.cnvkit.cnr -o proband_gene_cnv.tsv
  %(prog)s input.cnr output.tsv --del-threshold -0.4 --dup-threshold 0.3
  %(prog)s input.cnr output.tsv --alpha 0.05 --min-bins 2

Algorithm:
  1. CBS (Circular Binary Segmentation) identifies copy number change points
  2. Segments are mapped to genes based on genomic overlap
  3. Copy number (CN) is calculated: CN = 2 * 2^log2
  4. CNV call (DEL/DUP/Normal) is based on configurable thresholds
  5. Confidence level reflects statistical significance and data quality

Thresholds:
  --del-threshold: log2 < this value is called DEL (default: -0.4)
                   Heterozygous deletion CN=1 has log2 ≈ -0.5
  --dup-threshold: log2 > this value is called DUP (default: 0.3)
                   Heterozygous duplication CN=3 has log2 ≈ 0.58
        """
    )
    parser.add_argument("-i", "--input", required=True, help="Input CNVkit .cnr file")
    parser.add_argument("-o", "--output", required=True, help="Output gene-level CNV file")
    parser.add_argument("--segments", help="Output CBS segments file (optional)")
    parser.add_argument("--alpha", type=float, default=0.01,
                        help="Significance threshold for CBS (default: 0.01)")
    parser.add_argument("--min-bins", type=int, default=3,
                        help="Minimum bins per segment (default: 3)")
    parser.add_argument("--del-threshold", type=float, default=-0.4,
                        help="Log2 threshold for deletion call (default: -0.4)")
    parser.add_argument("--dup-threshold", type=float, default=0.3,
                        help="Log2 threshold for duplication call (default: 0.3)")

    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Step 1: Read .cnr file
    print(f"Reading CNVkit .cnr file: {args.input}")
    bins = read_cnr_file(args.input)
    print(f"  Total target bins: {len(bins)}")

    # Step 2: Run CBS segmentation
    print(f"\nRunning CBS segmentation (alpha={args.alpha}, min_bins={args.min_bins})...")
    segments = run_cbs_segmentation(bins, args.alpha, args.min_bins)

    total_segments = sum(len(s) for s in segments.values())
    print(f"\nTotal segments: {total_segments}")

    # Step 3: Map to genes with thresholds
    print(f"\nMapping segments to genes...")
    print(f"  DEL threshold: {args.del_threshold}")
    print(f"  DUP threshold: {args.dup_threshold}")
    gene_cnvs = gene_level_cnv(bins, segments, args.del_threshold, args.dup_threshold)
    print(f"  Gene-level CNV calls: {len(gene_cnvs)}")

    # Step 4: Write output
    print(f"\nWriting output to: {args.output}")
    write_output(gene_cnvs, args.output)

    if args.segments:
        print(f"Writing segments to: {args.segments}")
        write_segments_file(segments, args.segments)

    # Summary
    print(f"\n" + "="*50)
    print("CNV Summary:")
    print("="*50)

    cnv_counts = Counter(c.cnv_call for c in gene_cnvs)
    conf_counts = Counter(c.confidence for c in gene_cnvs)

    print(f"  Normal: {cnv_counts.get('Normal', 0)}")
    print(f"  DEL (deletion): {cnv_counts.get('DEL', 0)}")
    print(f"  DUP (duplication): {cnv_counts.get('DUP', 0)}")
    print(f"\nConfidence levels:")
    print(f"  HIGH: {conf_counts.get('HIGH', 0)}")
    print(f"  MEDIUM: {conf_counts.get('MEDIUM', 0)}")
    print(f"  LOW: {conf_counts.get('LOW', 0)}")

    # CN distribution (rounded to nearest integer for summary)
    cn_counts = Counter(round(c.cn) for c in gene_cnvs)
    print(f"\nCopy number distribution:")
    for cn in sorted(cn_counts.keys()):
        print(f"  CN={cn}: {cn_counts[cn]}")

    # High-confidence CNVs
    hc_del = [c for c in gene_cnvs if c.cnv_call == "DEL" and c.confidence == "HIGH"]
    hc_dup = [c for c in gene_cnvs if c.cnv_call == "DUP" and c.confidence == "HIGH"]

    if hc_del:
        print(f"\nHigh-confidence deletions ({len(hc_del)}):")
        for c in sorted(hc_del, key=lambda x: x.log2_mean)[:10]:
            print(f"  {c.gene} ({c.chromosome}:{c.start}-{c.end}): log2={c.log2_mean:.3f}, CN={c.cn:.2f}")

    if hc_dup:
        print(f"\nHigh-confidence duplications ({len(hc_dup)}):")
        for c in sorted(hc_dup, key=lambda x: -x.log2_mean)[:10]:
            print(f"  {c.gene} ({c.chromosome}:{c.start}-{c.end}): log2={c.log2_mean:.3f}, CN={c.cn:.2f}")

    print(f"\nDone!")


if __name__ == "__main__":
    main()