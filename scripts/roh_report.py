#!/usr/bin/env python3
"""
ROH (Runs of Homozygosity) Report Generator
Analyzes ROH regions for recessive disease genes overlap

Usage:
    python roh_report.py -i proband.ROH.txt -o roh_report.txt
"""

import argparse
import sys
import pandas as pd
from pathlib import Path


def load_recessive_genes(gencc_path):
    """Load recessive disease genes from GenCC submissions file"""
    try:
        df = pd.read_excel(gencc_path)
    except Exception as e:
        print(f"Error loading GenCC file: {e}")
        return set()

    # Filter for recessive inheritance modes
    recessive_modes = ['Autosomal recessive', 'X-linked recessive']
    recessive_df = df[df['moi_title'].isin(recessive_modes)]

    # Get unique gene symbols
    recessive_genes = set(recessive_df['gene_symbol'].unique())

    print(f"Found {len(recessive_genes)} recessive disease genes from GenCC")
    return recessive_genes


def load_gene_regions(bed_path):
    """Load gene regions from Gencode BED file"""
    gene_regions = {}  # gene -> list of (chrom, start, end)

    try:
        with open(bed_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 4:
                    continue

                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])

                # Parse gene info: gene|transcript|ensembl|exon|strand|cytoband
                info = parts[3]
                gene = info.split('|')[0] if '|' in info else info

                if gene not in gene_regions:
                    gene_regions[gene] = []
                gene_regions[gene].append((chrom, start, end))

    except Exception as e:
        print(f"Error loading BED file: {e}")
        return {}

    print(f"Found {len(gene_regions)} genes in BED file")
    return gene_regions


def get_gene_span(gene_regions, gene):
    """Get the full span of a gene (min start to max end)"""
    if gene not in gene_regions:
        return None

    regions = gene_regions[gene]
    if not regions:
        return None

    # Get the chromosome (should be consistent)
    chrom = regions[0][0]

    # Get min start and max end
    min_start = min(r[1] for r in regions)
    max_end = max(r[2] for r in regions)

    return (chrom, min_start, max_end)


def regions_overlap(chrom1, start1, end1, chrom2, start2, end2):
    """Check if two regions overlap"""
    if chrom1 != chrom2:
        return False

    # Overlap if: start1 < end2 and start2 < end1
    return start1 < end2 and start2 < end1


def find_recessive_genes_in_roh(roh_chrom, roh_start, roh_end,
                                 recessive_genes, gene_regions):
    """Find recessive disease genes that overlap with ROH region"""
    overlapping_genes = []

    for gene in recessive_genes:
        gene_span = get_gene_span(gene_regions, gene)
        if gene_span is None:
            continue

        gene_chrom, gene_start, gene_end = gene_span

        # Convert chromosome format (ROH uses numeric, BED uses numeric)
        # Handle chr prefix if present
        if gene_chrom.startswith('chr'):
            gene_chrom = gene_chrom[3:]

        if regions_overlap(roh_chrom, roh_start, roh_end,
                          gene_chrom, gene_start, gene_end):
            overlapping_genes.append(gene)

    return sorted(overlapping_genes)


def parse_roh_file(roh_path):
    """Parse ROH file and return list of records"""
    records = []

    with open(roh_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#Chr'):
                continue

            parts = line.split('\t')
            if len(parts) < 6:
                continue

            record = {
                'chrom': parts[0],
                'begin': int(parts[1]),
                'end': int(parts[2]),
                'size_mb': float(parts[3]),
                'nb_variants': int(parts[4]),
                'percent_homo': float(parts[5])
            }
            records.append(record)

    return records


def generate_report(roh_records, recessive_genes, gene_regions, output_path):
    """Generate ROH report with recessive genes overlap"""

    headers = [
        'Chr',
        'Begin',
        'End',
        'Size(Mb)',
        'Nb_variants',
        'Percentage_homozygosity',
        'Recessive_Genes',
        'Gene_Count'
    ]

    total_genes_found = set()

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(headers) + '\n')

        for record in roh_records:
            # Find overlapping recessive genes
            overlapping = find_recessive_genes_in_roh(
                record['chrom'],
                record['begin'],
                record['end'],
                recessive_genes,
                gene_regions
            )

            # Format genes output
            genes_str = ','.join(overlapping) if overlapping else '.'
            gene_count = len(overlapping)

            # Track all genes found
            total_genes_found.update(overlapping)

            row = [
                record['chrom'],
                record['begin'],
                record['end'],
                record['size_mb'],
                record['nb_variants'],
                record['percent_homo'],
                genes_str,
                gene_count
            ]

            f.write('\t'.join(str(x) for x in row) + '\n')

    return total_genes_found


def generate_statistics(roh_records, total_genes_found):
    """Generate summary statistics"""
    stats = {
        'total_roh_regions': len(roh_records),
        'total_size_mb': sum(r['size_mb'] for r in roh_records),
        'regions_with_genes': 0,
        'total_recessive_genes': len(total_genes_found),
        'largest_region': None,
        'genes_by_region': {}
    }

    for record in roh_records:
        if record['size_mb'] > (stats['largest_region']['size_mb'] if stats['largest_region'] else 0):
            stats['largest_region'] = record

    return stats


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='ROH Report Generator - Analyze recessive disease genes in ROH regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example usage:
  python roh_report.py -i proband.ROH.txt -o roh_report.txt
  python roh_report.py -i proband.ROH.txt -o roh_report.txt -g gencc.xlsx -b gencode.bed
        '''
    )

    parser.add_argument('-i', '--input',
                        required=True,
                        help='Input ROH file')

    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output report file path')

    parser.add_argument('-g', '--gencc',
                        default=None,
                        help='GenCC submissions Excel file (default: assets/gencc-submissions.xlsx)')

    parser.add_argument('-b', '--bed',
                        default=None,
                        help='Gencode BED file (default: assets/Gencode.GRCh37.cnvkit.target.bed)')

    args = parser.parse_args()

    # Set default paths
    base_path = Path(__file__).parent.parent
    gencc_path = Path(args.gencc) if args.gencc else base_path / 'assets' / 'gencc-submissions.xlsx'
    bed_path = Path(args.bed) if args.bed else base_path / 'assets' / 'Gencode.GRCh37.cnvkit.target.bed'

    roh_path = Path(args.input)
    output_path = Path(args.output)

    print(f"Input ROH: {roh_path}")
    print(f"Output report: {output_path}")
    print(f"GenCC file: {gencc_path}")
    print(f"BED file: {bed_path}")

    # Check files exist
    if not roh_path.exists():
        print(f"Error: ROH file not found: {roh_path}")
        sys.exit(1)

    if not gencc_path.exists():
        print(f"Error: GenCC file not found: {gencc_path}")
        sys.exit(1)

    if not bed_path.exists():
        print(f"Error: BED file not found: {bed_path}")
        sys.exit(1)

    # Load recessive genes
    print("\nLoading recessive disease genes...")
    recessive_genes = load_recessive_genes(gencc_path)

    # Load gene regions
    print("\nLoading gene regions from BED...")
    gene_regions = load_gene_regions(bed_path)

    # Parse ROH file
    print("\nParsing ROH file...")
    roh_records = parse_roh_file(roh_path)
    print(f"Found {len(roh_records)} ROH regions")

    # Generate report
    print("\nGenerating report...")
    total_genes = generate_report(roh_records, recessive_genes, gene_regions, output_path)
    print(f"Report written with {len(roh_records)} entries")

    # Summary statistics
    stats = generate_statistics(roh_records, total_genes)

    print("\n=== ROH Summary Statistics ===")
    print(f"Total ROH regions: {stats['total_roh_regions']}")
    print(f"Total ROH size: {stats['total_size_mb']:.2f} Mb")

    # Count regions with genes
    with open(output_path, 'r') as f:
        lines = f.readlines()[1:]  # Skip header
        regions_with_genes = sum(1 for line in lines
                                if line.strip().split('\t')[-2] != '.')
    print(f"Regions with recessive genes: {regions_with_genes}")
    print(f"Total recessive genes found: {stats['total_recessive_genes']}")

    if stats['largest_region']:
        print(f"Largest ROH region: {stats['largest_region']['chrom']}:"
              f"{stats['largest_region']['begin']}-{stats['largest_region']['end']} "
              f"({stats['largest_region']['size_mb']:.2f} Mb)")

    if total_genes:
        print(f"\nRecessive genes in ROH regions: {', '.join(sorted(total_genes))}")

    print(f"\nReport saved to: {output_path}")


if __name__ == '__main__':
    main()