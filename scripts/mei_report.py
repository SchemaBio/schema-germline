#!/usr/bin/env python3
"""
MEI (Mobile Element Insertion) Report Generator
Converts VEP-annotated MEI VCF to WES transposon analysis report format

Usage:
    python mei_report.py -i input.vcf -o output.txt
"""

import re
import sys
import json
import argparse
from pathlib import Path
from collections import defaultdict


# Path to transcripts configuration file
TRANSCRIPTS_FILE = Path(__file__).parent.parent / 'assets' / 'transcripts.json'

# Global transcripts lookup cache
_transcripts_cache = None


def load_transcripts():
    """Load transcripts configuration file"""
    global _transcripts_cache
    if _transcripts_cache is None:
        try:
            with open(TRANSCRIPTS_FILE, 'r', encoding='utf-8') as f:
                _transcripts_cache = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            _transcripts_cache = {}
    return _transcripts_cache


def strip_transcript_version(transcript_id):
    """Remove version number from transcript ID (e.g., NM_014110.5 -> NM_014110)"""
    if transcript_id and transcript_id != '.':
        # Split on dot and take first part
        return transcript_id.split('.')[0]
    return transcript_id


def get_preferred_transcript_id(gene):
    """Get preferred RefSeq transcript ID for a gene (without version)"""
    transcripts = load_transcripts()
    if gene in transcripts:
        refseq_id = transcripts[gene].get('refseq', '')
        return strip_transcript_version(refseq_id)
    return None


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


def parse_csq_field(csq_str, csq_header):
    """Parse VEP CSQ field into list of dictionaries"""
    if not csq_str or csq_str == '.':
        return []

    fields = csq_header.split('|')
    csq_entries = []

    for entry in csq_str.split(','):
        values = entry.split('|')
        csq_dict = {}
        for i, field in enumerate(fields):
            if i < len(values):
                csq_dict[field] = values[i] if values[i] else '.'
            else:
                csq_dict[field] = '.'
        csq_entries.append(csq_dict)

    return csq_entries


def get_impact_rank(impact):
    """Rank impact for sorting (HIGH=0, MODERATE=1, LOW=2, MODIFIER=3)"""
    impact_order = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
    return impact_order.get(impact, 4)


def format_gnomad_af(csq_dict):
    """Format GnomAD allele frequency from CSQ - returns AF value or '.' if empty"""
    af = csq_dict.get('GnomAD_AF_joint', '.')

    if af == '.' or af == '':
        return '.'

    try:
        # Handle multiple values separated by '&'
        af_val = float(af.split('&')[0]) if '&' in af else float(af)
        return f"{af_val:.6f}"
    except (ValueError, IndexError):
        return '.'


def get_best_csq_entry(csq_entries, gene=None):
    """Get the most significant CSQ entry with preferred transcript priority

    Selection logic:
    1. Among entries with highest IMPACT, prefer those matching transcripts.json
    2. If no preferred transcript found among highest impact entries,
       fall back to first entry with highest impact
    """
    if not csq_entries:
        return None

    # Sort by impact rank (HIGH=0, MODERATE=1, LOW=2, MODIFIER=3)
    sorted_entries = sorted(csq_entries, key=lambda x: get_impact_rank(x.get('IMPACT', 'MODIFIER')))

    # Get all entries with the highest impact
    highest_impact_rank = get_impact_rank(sorted_entries[0].get('IMPACT', 'MODIFIER'))
    highest_impact_entries = [e for e in sorted_entries
                              if get_impact_rank(e.get('IMPACT', 'MODIFIER')) == highest_impact_rank]

    # Get preferred transcript for the gene
    preferred_tx_id = None
    if gene:
        preferred_tx_id = get_preferred_transcript_id(gene)

    # Try to find an entry matching the preferred transcript
    if preferred_tx_id:
        for entry in highest_impact_entries:
            feature = entry.get('Feature', '.')
            if feature and feature != '.':
                feature_base = strip_transcript_version(feature)
                if feature_base == preferred_tx_id:
                    return entry

    # Fall back to first entry with highest impact
    return highest_impact_entries[0]


def extract_transcript_info(csq_dict):
    """Extract exon/intron information"""
    exon = csq_dict.get('EXON', '.')
    intron = csq_dict.get('INTRON', '.')
    distance = csq_dict.get('DISTANCE', '.')

    location = ''
    if exon and exon != '.':
        location = f"Exon {exon}"
    elif intron and intron != '.':
        location = f"Intron {intron}"
    elif distance and distance != '.':
        location = f"Distance: {distance}"

    return location


def format_consequence(consequence):
    """Format consequence for readability"""
    # Clean up consequence string
    consequence = consequence.replace('_', ' ')
    consequence = consequence.replace('&', ', ')
    return consequence


def parse_vcf(vcf_path):
    """Parse VCF file and return list of MEI records"""
    records = []
    csq_header = None

    with open(vcf_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()

            if line.startswith('##'):
                # Extract CSQ header format
                if 'CSQ' in line and 'Format:' in line:
                    match = re.search(r'Format:\s*(.+)"', line)
                    if match:
                        csq_header = match.group(1)
                continue

            if line.startswith('#CHROM'):
                continue

            if not line:
                continue

            # Parse data line
            parts = line.split('\t')
            if len(parts) < 8:
                continue

            chrom = parts[0]
            pos = parts[1]
            var_id = parts[2]
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]
            filter_status = parts[6]
            info_str = parts[7]

            info = parse_info_field(info_str)

            # Parse CSQ
            csq_str = info.get('CSQ', '')
            csq_entries = parse_csq_field(csq_str, csq_header) if csq_header else []

            record = {
                'chrom': chrom,
                'pos': pos,
                'var_id': var_id,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'filter': filter_status,
                'info': info,
                'csq_entries': csq_entries
            }
            records.append(record)

    return records, csq_header


def generate_report(records, output_path):
    """Generate MEI report in tab-separated format"""

    # Define column headers
    headers = [
        'Chromosome',
        'Position',
        'MEI_ID',
        'TE_Type',
        'TE_Family',
        'Direction',
        'Confidence',
        'Support_Reads',
        'Avg_SoftClip_Length',
        'Gene',
        'Transcript',
        'Location',
        'Consequence',
        'Impact',
        'Cytoband',
        'ClinVar_Sig',
        'ClinVar_DN',
        'ClinVar_Star',
        'GnomAD_AF',
        'HGNC_ID',
        'Filter'
    ]

    # Sort records by confidence and impact
    def sort_key(r):
        conf_order = {'High': 0, 'Medium': 1, 'Low': 2}
        conf = conf_order.get(r['info'].get('CONFIDENCE', 'Low'), 2)

        # Get first gene symbol from CSQ entries for transcript selection
        first_gene = None
        if r['csq_entries']:
            first_gene = r['csq_entries'][0].get('SYMBOL', None)

        best_csq = get_best_csq_entry(r['csq_entries'], gene=first_gene)
        impact = get_impact_rank(best_csq.get('IMPACT', 'MODIFIER') if best_csq else 'MODIFIER')

        # Sort by chromosome and position as secondary keys
        try:
            chrom_num = int(r['chrom'])
        except ValueError:
            chrom_num = 100  # X, Y, MT etc
        pos = int(r['pos'])

        return (conf, impact, chrom_num, pos)

    sorted_records = sorted(records, key=sort_key)

    # Write report
    with open(output_path, 'w', encoding='utf-8') as f:
        # Write header
        f.write('\t'.join(headers) + '\n')

        # Write data rows
        for record in sorted_records:
            info = record['info']

            # Get first gene symbol from CSQ entries for transcript selection
            first_gene = None
            if record['csq_entries']:
                first_gene = record['csq_entries'][0].get('SYMBOL', None)

            best_csq = get_best_csq_entry(record['csq_entries'], gene=first_gene)

            # Extract TE type from ALT or MINAME
            te_type = info.get('MINAME', 'UNKNOWN')
            if 'ALU' in record['alt']:
                te_type = 'Alu'
            elif 'LINE1' in record['alt']:
                te_type = 'L1'
            elif 'SVA' in record['alt']:
                te_type = 'SVA'
            elif 'HERV' in record['alt']:
                te_type = 'HERV'

            # Format direction
            direction = info.get('DIR', '.')
            if direction == '5_prime':
                direction = '5\''
            elif direction == '3_prime':
                direction = '3\''

            # Format confidence
            confidence = info.get('CONFIDENCE', '.')

            # Support reads
            support = info.get('SUPPORT', '.')

            # Avg softclip length
            avg_sclen = info.get('AVGSCLEN', '.')

            # TE family
            te_family = info.get('MEFAMILY', '.')

            # Gene and transcript info from CSQ
            if best_csq:
                gene = best_csq.get('SYMBOL', '.')
                transcript = best_csq.get('Feature', '.')
                location = extract_transcript_info(best_csq)
                consequence = format_consequence(best_csq.get('Consequence', '.'))
                impact = best_csq.get('IMPACT', '.')
                cytoband = best_csq.get('cytoBand', '.')
                clinvar_sig = best_csq.get('CLNSIG', '.') if best_csq.get('CLNSIG', '.') else '.'
                clinvar_dn = best_csq.get('CLNDN', '.') if best_csq.get('CLNDN', '.') else '.'
                clinvar_star = best_csq.get('CLNSTAR', '.') if best_csq.get('CLNSTAR', '.') else '.'
                hgnc_id = best_csq.get('HGNC_ID', '.')
                gnomad_af = format_gnomad_af(best_csq)
            else:
                gene = '.'
                transcript = '.'
                location = 'Intergenic'
                consequence = 'mobile_element_insertion'
                impact = 'MODIFIER'
                cytoband = '.'
                clinvar_sig = '.'
                clinvar_dn = '.'
                clinvar_star = '.'
                hgnc_id = '.'
                gnomad_af = '.'

            row = [
                record['chrom'],
                record['pos'],
                record['var_id'],
                te_type,
                te_family,
                direction,
                confidence,
                support,
                avg_sclen,
                gene,
                transcript,
                location,
                consequence,
                impact,
                cytoband,
                clinvar_sig,
                clinvar_dn,
                clinvar_star,
                gnomad_af,
                hgnc_id,
                record['filter']
            ]

            f.write('\t'.join(row) + '\n')

    return len(sorted_records)


def generate_summary_statistics(records):
    """Generate summary statistics for MEI report"""
    stats = {
        'total': len(records),
        'pass': 0,
        'low_confidence': 0,
        'by_type': defaultdict(int),
        'by_impact': defaultdict(int),
        'high_impact_genes': []
    }

    for record in records:
        info = record['info']

        # Count by filter
        if record['filter'] == 'PASS':
            stats['pass'] += 1
        else:
            stats['low_confidence'] += 1

        # Count by TE type
        te_type = info.get('MINAME', 'UNKNOWN')
        stats['by_type'][te_type] += 1

        # Get first gene symbol for transcript selection
        first_gene = None
        if record['csq_entries']:
            first_gene = record['csq_entries'][0].get('SYMBOL', None)

        # Count by impact
        best_csq = get_best_csq_entry(record['csq_entries'], gene=first_gene)
        if best_csq:
            impact = best_csq.get('IMPACT', 'MODIFIER')
            stats['by_impact'][impact] += 1

            if impact == 'HIGH' and best_csq.get('SYMBOL', '.') != '.':
                stats['high_impact_genes'].append(best_csq.get('SYMBOL'))

    return stats


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='MEI (Mobile Element Insertion) Report Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example usage:
  python mei_report.py -i proband.mei.vep.vcf -o mei_report.txt
  python mei_report.py -i proband.mei.vep.vcf -o mei_report.txt -t transcripts.json
        '''
    )

    parser.add_argument('-i', '--input',
                        required=True,
                        help='Input VEP-annotated MEI VCF file')

    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output report file path')

    parser.add_argument('-t', '--transcripts',
                        default=None,
                        help='Transcripts configuration file (default: assets/transcripts.json)')

    args = parser.parse_args()

    # Set paths
    vcf_path = Path(args.input)
    output_path = Path(args.output)

    # Handle transcripts file
    if args.transcripts:
        global TRANSCRIPTS_FILE
        TRANSCRIPTS_FILE = Path(args.transcripts)

    print(f"Input VCF: {vcf_path}")
    print(f"Output report: {output_path}")
    if args.transcripts:
        print(f"Transcripts file: {TRANSCRIPTS_FILE}")

    if not vcf_path.exists():
        print(f"Error: VCF file not found: {vcf_path}")
        sys.exit(1)

    # Parse VCF
    print("Parsing VCF file...")
    records, csq_header = parse_vcf(vcf_path)
    print(f"Found {len(records)} MEI records")

    if not csq_header:
        print("Warning: CSQ header not found, annotation data may be incomplete")

    # Generate report
    print("Generating report...")
    count = generate_report(records, output_path)
    print(f"Report written with {count} entries")

    # Generate summary
    stats = generate_summary_statistics(records)
    print("\n=== MEI Summary Statistics ===")
    print(f"Total MEI events: {stats['total']}")
    print(f"PASS (High/Medium confidence): {stats['pass']}")
    print(f"Filtered (Low confidence): {stats['low_confidence']}")
    print("\nBy TE Type:")
    for te_type, count in sorted(stats['by_type'].items()):
        print(f"  {te_type}: {count}")
    print("\nBy Impact:")
    for impact, count in sorted(stats['by_impact'].items(), key=lambda x: get_impact_rank(x[0])):
        print(f"  {impact}: {count}")
    if stats['high_impact_genes']:
        print(f"\nHIGH impact genes: {', '.join(set(stats['high_impact_genes']))}")

    print(f"\nReport saved to: {output_path}")


if __name__ == '__main__':
    main()