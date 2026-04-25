#!/usr/bin/env python3
"""
STR (Short Tandem Repeat) Report Generator
Converts STR VCF to WES dynamic mutation analysis report format

Usage:
    python str_report.py -i input.str.anno.vcf -o output.txt
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict


def parse_info_field(info_str):
    """Parse VCF INFO field into dictionary"""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict


def parse_format_field(format_str, sample_str):
    """Parse FORMAT and sample fields into dictionary"""
    format_fields = format_str.split(':')
    sample_values = sample_str.split(':')

    format_dict = {}
    for i, field in enumerate(format_fields):
        if i < len(sample_values):
            format_dict[field] = sample_values[i]
        else:
            format_dict[field] = '.'
    return format_dict


def parse_gt_repcn(gt, repcn):
    """Parse genotype and repeat counts into allele information"""
    gt_parts = gt.split('/')
    repcn_parts = repcn.split('/') if repcn and repcn != '.' else ['.']

    alleles = []
    for i in range(len(gt_parts)):
        allele_gt = gt_parts[i] if i < len(gt_parts) else '.'
        allele_cn = repcn_parts[i] if i < len(repcn_parts) else '.'
        alleles.append({
            'allele': allele_gt,
            'repeat_count': allele_cn
        })
    return alleles


def get_str_status_display(status, filter_status):
    """Convert STR status to display format"""
    if not status or status == '.':
        # If no status and filter is LowDepth, mark as QC failed
        if filter_status != 'PASS':
            return 'QC_Failed'
        return 'QC_Failed'

    status_map = {
        'normal': 'Normal',
        'pre_mutation': 'Pre-mutation',
        'full_mutation': 'Full mutation'
    }
    return status_map.get(status, status)


def evaluate_pathogenicity(alleles, normal_max, pathologic_min, inheritance, filter_status):
    """Evaluate pathogenicity based on repeat counts and thresholds"""

    # If QC failed, return QC failed
    if filter_status != 'PASS':
        return 'QC_Failed'

    # If no thresholds available
    if not normal_max or not pathologic_min or normal_max == '.' or pathologic_min == '.':
        return 'No threshold data'

    try:
        normal_max_val = int(normal_max)
        pathologic_min_val = int(pathologic_min)
    except ValueError:
        return 'No threshold data'

    # Check each allele
    pathologic_alleles = []
    pre_mutation_alleles = []

    for i, allele in enumerate(alleles):
        cn = allele['repeat_count']
        if cn == '.' or not cn:
            return 'Unknown (no data)'

        try:
            cn_val = int(cn)
        except ValueError:
            return 'Unknown (invalid data)'

        if cn_val >= pathologic_min_val:
            pathologic_alleles.append(i + 1)  # Allele 1 or 2
        elif cn_val > normal_max_val:
            pre_mutation_alleles.append(i + 1)

    # Determine final interpretation
    if pathologic_alleles:
        if inheritance == 'AD':
            # Autosomal dominant: one pathologic allele is enough
            return 'Pathogenic'
        elif inheritance == 'AR':
            # Autosomal recessive: need two pathologic alleles
            if len(pathologic_alleles) >= 2:
                return 'Pathogenic'
            else:
                return 'Likely pathogenic (carrier)'
        elif inheritance in ['XR', 'XD']:
            # X-linked: depends on sex, but generally pathogenic if present
            return 'Pathogenic'
        return 'Pathogenic'

    if pre_mutation_alleles:
        if len(pre_mutation_alleles) >= 2:
            return 'Pre-mutation (both alleles)'
        return f'Pre-mutation (allele {pre_mutation_alleles[0]})'

    return 'Benign'


def get_inheritance_display(mode):
    """Convert inheritance mode to display format"""
    mode_map = {
        'AD': 'Autosomal Dominant',
        'AR': 'Autosomal Recessive',
        'XR': 'X-linked Recessive',
        'XD': 'X-linked Dominant'
    }
    return mode_map.get(mode, mode)


def format_repcn_display(alleles, ref_cn, normal_max, pathologic_min):
    """Format repeat count display with interpretation"""
    display_parts = []
    for allele in alleles:
        cn = allele['repeat_count']
        if cn == '.':
            display_parts.append('NA')
        else:
            try:
                cn_val = int(cn)
                # Determine status based on thresholds
                if normal_max and pathologic_min:
                    try:
                        normal_max_val = int(normal_max)
                        pathologic_min_val = int(pathologic_min)
                        if cn_val <= normal_max_val:
                            status = 'N'  # Normal
                        elif cn_val >= pathologic_min_val:
                            status = 'P'  # Pathologic
                        else:
                            status = 'PM'  # Pre-mutation
                        display_parts.append(f"{cn_val}({status})")
                    except ValueError:
                        display_parts.append(str(cn_val))
                else:
                    display_parts.append(str(cn_val))
            except ValueError:
                display_parts.append(cn)
    return '/'.join(display_parts)


def parse_vcf(vcf_path):
    """Parse STR VCF file and return list of records"""
    records = []

    with open(vcf_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()

            if line.startswith('##'):
                continue

            if line.startswith('#CHROM'):
                # Header line - get sample name
                parts = line.split('\t')
                sample_name = parts[-1] if len(parts) > 9 else 'sample'
                continue

            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 10:
                continue

            chrom = parts[0]
            pos = parts[1]
            var_id = parts[2]
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]
            filter_status = parts[6]
            info_str = parts[7]
            format_str = parts[8]
            sample_str = parts[9]

            info = parse_info_field(info_str)
            format_data = parse_format_field(format_str, sample_str)

            # Parse genotype and repeat counts
            gt = format_data.get('GT', '.')
            repcn = format_data.get('REPCN', '.')
            alleles = parse_gt_repcn(gt, repcn)

            record = {
                'chrom': chrom,
                'pos': pos,
                'var_id': var_id,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'filter': filter_status,
                'info': info,
                'format': format_data,
                'alleles': alleles
            }
            records.append(record)

    return records


def generate_report(records, output_path):
    """Generate STR report in tab-separated format"""

    # Define column headers
    headers = [
        'Chromosome',
        'Position',
        'Gene',
        'Repeat_Unit',
        'Ref_Repeats',
        'Allele1_Repeats',
        'Allele2_Repeats',
        'Repeat_Display',
        'STR_Status',
        'Pathogenicity',
        'Normal_Max',
        'Pathologic_Min',
        'Disease',
        'Inheritance',
        'HGNC_ID',
        'Depth',
        'Spanning_Reads',
        'Flanking_Reads',
        'InRepeat_Reads',
        'SweGen_Mean',
        'SweGen_Std',
        'Source',
        'Filter'
    ]

    # Sort records by filter status and STR status
    def sort_key(r):
        # Priority: PASS over LowDepth
        filter_order = {'PASS': 0, 'LowDepth': 1}
        filter_val = filter_order.get(r['filter'], 2)

        # Secondary: by gene name
        gene = r['info'].get('REPID', '')

        return (filter_val, gene)

    sorted_records = sorted(records, key=sort_key)

    # Write report
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(headers) + '\n')

        for record in sorted_records:
            info = record['info']
            format_data = record['format']
            alleles = record['alleles']

            # Gene/Repeat ID
            gene = info.get('REPID', '.')

            # Repeat unit (use display format if available)
            repeat_unit = info.get('DisplayRU', info.get('RU', '.'))

            # Reference repeats
            ref_repeats = info.get('REF', '.')

            # Allele repeat counts
            allele1_cn = alleles[0]['repeat_count'] if alleles else '.'
            allele2_cn = alleles[1]['repeat_count'] if len(alleles) > 1 else '.'

            # Format display with interpretation
            repeat_display = format_repcn_display(
                alleles,
                ref_repeats,
                info.get('STR_NORMAL_MAX', ''),
                info.get('STR_PATHOLOGIC_MIN', '')
            )

            # STR status
            str_status = get_str_status_display(info.get('STR_STATUS', '.'), record['filter'])

            # Pathogenicity evaluation
            pathogenicity = evaluate_pathogenicity(
                alleles,
                info.get('STR_NORMAL_MAX', '.'),
                info.get('STR_PATHOLOGIC_MIN', '.'),
                info.get('InheritanceMode', '.'),
                record['filter']
            )

            # Thresholds
            normal_max = info.get('STR_NORMAL_MAX', '.')
            pathologic_min = info.get('STR_PATHOLOGIC_MIN', '.')

            # Disease
            disease = info.get('Disease', '.')

            # Inheritance mode
            inheritance = get_inheritance_display(info.get('InheritanceMode', '.'))

            # HGNC ID
            hgnc_id = info.get('HGNCId', '.')

            # Depth/coverage
            depth = format_data.get('LC', '.')

            # Supporting reads
            adsp = format_data.get('ADSP', '.')
            adfl = format_data.get('ADFL', '.')
            adir = format_data.get('ADIR', '.')

            # Format reads display (allele1/allele2)
            spanning_reads = adsp if adsp != '.' else '.'
            flanking_reads = adfl if adfl != '.' else '.'
            inrepeat_reads = adir if adir != '.' else '.'

            # Population data
            swegen_mean = info.get('SweGenMean', '.')
            swegen_std = info.get('SweGenStd', '.')

            # Source
            source = info.get('SourceDisplay', info.get('Source', '.'))

            row = [
                record['chrom'],
                record['pos'],
                gene,
                repeat_unit,
                ref_repeats,
                allele1_cn,
                allele2_cn,
                repeat_display,
                str_status,
                pathogenicity,
                normal_max,
                pathologic_min,
                disease,
                inheritance,
                hgnc_id,
                depth,
                spanning_reads,
                flanking_reads,
                inrepeat_reads,
                swegen_mean,
                swegen_std,
                source,
                record['filter']
            ]

            f.write('\t'.join(row) + '\n')

    return len(sorted_records)


def generate_summary_statistics(records):
    """Generate summary statistics for STR report"""
    stats = {
        'total': len(records),
        'pass': 0,
        'low_depth': 0,
        'by_status': defaultdict(int),
        'by_inheritance': defaultdict(int),
        'diseases': set()
    }

    for record in records:
        info = record['info']

        # Count by filter
        if record['filter'] == 'PASS':
            stats['pass'] += 1
        else:
            stats['low_depth'] += 1

        # Count by STR status
        status = info.get('STR_STATUS', 'unknown')
        if status:
            stats['by_status'][status] += 1

        # Count by inheritance mode
        inheritance = info.get('InheritanceMode', '.')
        if inheritance and inheritance != '.':
            stats['by_inheritance'][inheritance] += 1

        # Collect diseases
        disease = info.get('Disease', '.')
        if disease and disease != '.':
            stats['diseases'].add(disease)

    return stats


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='STR (Short Tandem Repeat) Report Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example usage:
  python str_report.py -i proband.str.anno.vcf -o str_report.txt

Report columns:
  - Repeat_Display: Shows repeat counts with interpretation
    N = Normal, PM = Pre-mutation, P = Pathologic
    Example: 28(N)/29(N) means both alleles are in normal range
        '''
    )

    parser.add_argument('-i', '--input',
                        required=True,
                        help='Input STR VCF file')

    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output report file path')

    args = parser.parse_args()

    vcf_path = Path(args.input)
    output_path = Path(args.output)

    print(f"Input VCF: {vcf_path}")
    print(f"Output report: {output_path}")

    if not vcf_path.exists():
        print(f"Error: VCF file not found: {vcf_path}")
        sys.exit(1)

    # Parse VCF
    print("Parsing VCF file...")
    records = parse_vcf(vcf_path)
    print(f"Found {len(records)} STR records")

    # Generate report
    print("Generating report...")
    count = generate_report(records, output_path)
    print(f"Report written with {count} entries")

    # Generate summary
    stats = generate_summary_statistics(records)
    print("\n=== STR Summary Statistics ===")
    print(f"Total STR loci: {stats['total']}")
    print(f"PASS: {stats['pass']}")
    print(f"LowDepth (filtered): {stats['low_depth']}")
    print("\nBy STR Status:")
    for status, count in sorted(stats['by_status'].items()):
        print(f"  {status}: {count}")
    print("\nBy Inheritance Mode:")
    for mode, count in sorted(stats['by_inheritance'].items()):
        print(f"  {mode}: {count}")
    if stats['diseases']:
        print(f"\nDiseases tested: {', '.join(sorted(stats['diseases']))}")

    print(f"\nReport saved to: {output_path}")


if __name__ == '__main__':
    main()