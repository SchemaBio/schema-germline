#!/usr/bin/env python3
"""
UPD (Uniparental Disomy) Detection Script
Detects potential UPD regions by analyzing ROH overlap in trio (proband, father, mother)

Usage:
    python upd_detection.py -p proband_report.txt -f father_report.txt -m mother_report.txt -o upd_report.txt
"""

import argparse
import sys
import pandas as pd
from pathlib import Path


def parse_roh_report(file_path):
    """Parse ROH report file and return list of records

    Expected columns: Chr, Begin, End, Size(Mb), Nb_variants,
                      Percentage_homozygosity, Recessive_Genes, Gene_Count
    """
    records = []

    with open(file_path, 'r', encoding='utf-8') as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')

            # Detect header line
            if parts[0] == 'Chr':
                header = parts
                continue

            if len(parts) < 6:
                continue

            record = {
                'chrom': parts[0],
                'begin': int(parts[1]),
                'end': int(parts[2]),
                'size_mb': float(parts[3]),
                'nb_variants': int(parts[4]) if parts[4] != '.' else 0,
                'percent_homo': float(parts[5]) if parts[5] != '.' else 0,
                'recessive_genes': parts[6] if len(parts) > 6 and parts[6] != '.' else '',
                'gene_count': int(parts[7]) if len(parts) > 7 and parts[7] != '.' else 0
            }
            records.append(record)

    return records


def regions_overlap(chrom1, start1, end1, chrom2, start2, end2):
    """Check if two regions overlap"""
    if chrom1 != chrom2:
        return False

    # Normalize chromosome format
    if chrom1.startswith('chr'):
        chrom1 = chrom1[3:]
    if chrom2.startswith('chr'):
        chrom2 = chrom2[3:]

    if chrom1 != chrom2:
        return False

    # Overlap if: start1 < end2 and start2 < end1
    return start1 < end2 and start2 < end1


def find_overlapping_regions(proband_region, parent_regions):
    """Find all parent regions that overlap with proband region"""
    overlapping = []

    for region in parent_regions:
        if regions_overlap(
            proband_region['chrom'], proband_region['begin'], proband_region['end'],
            region['chrom'], region['begin'], region['end']
        ):
            overlapping.append(region)

    return overlapping


def calculate_overlap_size(chrom1, start1, end1, chrom2, start2, end2):
    """Calculate the size of overlap between two regions"""
    if chrom1 != chrom2:
        return 0

    # Normalize chromosome format
    if chrom1.startswith('chr'):
        chrom1 = chrom1[3:]
    if chrom2.startswith('chr'):
        chrom2 = chrom2[3:]

    if chrom1 != chrom2:
        return 0

    # Calculate overlap
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    if overlap_start < overlap_end:
        return overlap_end - overlap_start
    return 0


def detect_upd(proband_records, father_records, mother_records, min_size_mb=1.0):
    """
    Detect potential UPD regions

    UPD detection logic:
    - Proband has ROH in a region
    - Neither parent has ROH in the same region (both parents are heterozygous)
    - OR ROH only overlaps with one parent (suggesting isodisomy from that parent)

    Returns list of UPD candidate regions
    """
    upd_candidates = []

    for proband_region in proband_records:
        # Skip small regions
        if proband_region['size_mb'] < min_size_mb:
            continue

        chrom = proband_region['chrom']
        begin = proband_region['begin']
        end = proband_region['end']

        # Find overlapping regions in each parent
        father_overlaps = find_overlapping_regions(proband_region, father_records)
        mother_overlaps = find_overlapping_regions(proband_region, mother_records)

        # Determine UPD type
        father_has_roh = len(father_overlaps) > 0
        mother_has_roh = len(mother_overlaps) > 0

        upd_type = None
        confidence = 'Low'

        if not father_has_roh and not mother_has_roh:
            # Neither parent has ROH in this region
            # This is the strongest indicator of UPD
            # Could be either paternal or maternal isodisomy
            # Need additional SNP data to determine which parent
            upd_type = 'UPD (Unknown origin)'
            confidence = 'High'

        elif father_has_roh and not mother_has_roh:
            # Only father has ROH overlapping with proband
            # Could be paternal isodisomy (UPD-pat)
            # But this could also be normal inheritance where proband inherited
            # the same ROH region from father
            upd_type = 'Potential Paternal UPD (Isodisomy)'
            confidence = 'Medium'

        elif not father_has_roh and mother_has_roh:
            # Only mother has ROH overlapping with proband
            # Could be maternal isodisomy (UPD-mat)
            upd_type = 'Potential Maternal UPD (Isodisomy)'
            confidence = 'Medium'

        else:
            # Both parents have ROH overlapping with proband
            # This is less likely to be UPD, could be normal inheritance
            # or consanguinity
            upd_type = 'Unlikely UPD (Both parents have ROH)'
            confidence = 'Low'

        # Calculate overlap details
        father_overlap_size = 0
        mother_overlap_size = 0

        if father_has_roh:
            for fo in father_overlaps:
                father_overlap_size += calculate_overlap_size(
                    chrom, begin, end,
                    fo['chrom'], fo['begin'], fo['end']
                )

        if mother_has_roh:
            for mo in mother_overlaps:
                mother_overlap_size += calculate_overlap_size(
                    chrom, begin, end,
                    mo['chrom'], mo['begin'], mo['end']
                )

        candidate = {
            'chrom': chrom,
            'begin': begin,
            'end': end,
            'size_mb': proband_region['size_mb'],
            'percent_homo': proband_region['percent_homo'],
            'recessive_genes': proband_region['recessive_genes'],
            'gene_count': proband_region['gene_count'],
            'father_overlaps': len(father_overlaps),
            'mother_overlaps': len(mother_overlaps),
            'father_overlap_bp': father_overlap_size,
            'mother_overlap_bp': mother_overlap_size,
            'upd_type': upd_type,
            'confidence': confidence
        }

        upd_candidates.append(candidate)

    return upd_candidates


def generate_upd_report(upd_candidates, output_path, min_confidence='Low'):
    """Generate UPD detection report"""

    confidence_order = {'High': 0, 'Medium': 1, 'Low': 2}

    # Filter by confidence level
    min_level = confidence_order.get(min_confidence, 2)
    filtered_candidates = [
        c for c in upd_candidates
        if confidence_order.get(c['confidence'], 2) <= min_level
    ]

    # Sort by confidence (high first), then by size (large first)
    filtered_candidates.sort(
        key=lambda x: (confidence_order.get(x['confidence'], 2), -x['size_mb'])
    )

    headers = [
        'Chr',
        'Begin',
        'End',
        'Size(Mb)',
        'Homozygosity(%)',
        'Father_ROH_Overlap',
        'Mother_ROH_Overlap',
        'Father_Overlap(bp)',
        'Mother_Overlap(bp)',
        'Recessive_Genes',
        'Gene_Count',
        'UPD_Type',
        'Confidence'
    ]

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(headers) + '\n')

        for candidate in filtered_candidates:
            row = [
                candidate['chrom'],
                candidate['begin'],
                candidate['end'],
                f"{candidate['size_mb']:.2f}",
                f"{candidate['percent_homo']:.1f}",
                candidate['father_overlaps'],
                candidate['mother_overlaps'],
                candidate['father_overlap_bp'],
                candidate['mother_overlap_bp'],
                candidate['recessive_genes'] if candidate['recessive_genes'] else '.',
                candidate['gene_count'] if candidate['gene_count'] > 0 else '.',
                candidate['upd_type'],
                candidate['confidence']
            ]

            f.write('\t'.join(str(x) for x in row) + '\n')

    return filtered_candidates


def print_summary(upd_candidates):
    """Print summary statistics"""
    high_confidence = [c for c in upd_candidates if c['confidence'] == 'High']
    medium_confidence = [c for c in upd_candidates if c['confidence'] == 'Medium']

    print("\n" + "=" * 50)
    print("UPD Detection Summary")
    print("=" * 50)

    print(f"\nTotal proband ROH regions analyzed: {len(upd_candidates)}")
    print(f"High confidence UPD candidates: {len(high_confidence)}")
    print(f"Medium confidence UPD candidates: {len(medium_confidence)}")

    if high_confidence:
        print("\n--- High Confidence UPD Regions ---")
        for c in high_confidence:
            print(f"  {c['chrom']}:{c['begin']}-{c['end']} ({c['size_mb']:.2f} Mb)")
            if c['recessive_genes']:
                print(f"    Genes: {c['recessive_genes']}")

    if medium_confidence:
        print("\n--- Medium Confidence UPD Regions ---")
        for c in medium_confidence:
            print(f"  {c['chrom']}:{c['begin']}-{c['end']} ({c['size_mb']:.2f} Mb) - {c['upd_type']}")
            if c['recessive_genes']:
                print(f"    Genes: {c['recessive_genes']}")

    # Summary by UPD type
    print("\n--- Classification Summary ---")
    type_counts = {}
    for c in upd_candidates:
        t = c['upd_type']
        type_counts[t] = type_counts.get(t, 0) + 1

    for t, count in sorted(type_counts.items()):
        print(f"  {t}: {count}")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='UPD Detection - Analyze ROH patterns in trio for uniparental disomy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
UPD Detection Logic:
  - High confidence: Proband has ROH, but neither parent has ROH in the same region
  - Medium confidence: Proband has ROH overlapping with only one parent

Example usage:
  python upd_detection.py -p proband_report.txt -f father_report.txt -m mother_report.txt -o upd_report.txt
  python upd_detection.py -p proband_report.txt -f father_report.txt -m mother_report.txt -o upd_report.txt --min-size 2.0
        '''
    )

    parser.add_argument('-p', '--proband',
                        required=True,
                        help='Proband ROH report file')

    parser.add_argument('-f', '--father',
                        required=True,
                        help='Father ROH report file')

    parser.add_argument('-m', '--mother',
                        required=True,
                        help='Mother ROH report file')

    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output UPD report file path')

    parser.add_argument('--min-size',
                        type=float,
                        default=1.0,
                        help='Minimum ROH size in Mb to consider (default: 1.0)')

    parser.add_argument('--min-confidence',
                        choices=['High', 'Medium', 'Low'],
                        default='Low',
                        help='Minimum confidence level to report (default: Low)')

    args = parser.parse_args()

    proband_path = Path(args.proband)
    father_path = Path(args.father)
    mother_path = Path(args.mother)
    output_path = Path(args.output)

    print("UPD Detection Analysis")
    print("-" * 50)
    print(f"Proband: {proband_path}")
    print(f"Father: {father_path}")
    print(f"Mother: {mother_path}")
    print(f"Output: {output_path}")
    print(f"Min size: {args.min_size} Mb")
    print(f"Min confidence: {args.min_confidence}")

    # Check files exist
    for path, name in [(proband_path, 'Proband'), (father_path, 'Father'), (mother_path, 'Mother')]:
        if not path.exists():
            print(f"Error: {name} file not found: {path}")
            sys.exit(1)

    # Parse ROH reports
    print("\nParsing ROH reports...")
    proband_records = parse_roh_report(proband_path)
    father_records = parse_roh_report(father_path)
    mother_records = parse_roh_report(mother_path)

    print(f"  Proband: {len(proband_records)} ROH regions")
    print(f"  Father: {len(father_records)} ROH regions")
    print(f"  Mother: {len(mother_records)} ROH regions")

    # Detect UPD
    print("\nDetecting UPD regions...")
    upd_candidates = detect_upd(
        proband_records,
        father_records,
        mother_records,
        min_size_mb=args.min_size
    )

    # Generate report
    print("\nGenerating UPD report...")
    filtered = generate_upd_report(upd_candidates, output_path, args.min_confidence)
    print(f"Report written with {len(filtered)} entries")

    # Print summary
    print_summary(upd_candidates)

    print(f"\nReport saved to: {output_path}")


if __name__ == '__main__':
    main()