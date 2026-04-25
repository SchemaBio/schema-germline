#!/usr/bin/env python3
"""
VEP VCF Report Generator
Parses VEP-annotated VCF and generates a comprehensive variant report

Usage:
    python vep_report.py -i input.vcf -o output.txt
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
    """Format GnomAD allele frequency from CSQ"""
    af = csq_dict.get('GnomAD_AF_joint', '.')

    if af == '.' or af == '':
        return '.'

    try:
        af_val = float(af.split('&')[0]) if '&' in af else float(af)
        return f"{af_val:.6f}"
    except (ValueError, IndexError):
        return '.'


def format_gnomad_af_eas(csq_dict):
    """Format GnomAD East Asian allele frequency from CSQ"""
    af = csq_dict.get('GnomAD_AF_joint_eas', '.')

    if af == '.' or af == '':
        return '.'

    try:
        af_val = float(af.split('&')[0]) if '&' in af else float(af)
        return f"{af_val:.6f}"
    except (ValueError, IndexError):
        return '.'


def get_clinvar_field(csq_dict, field_name):
    """Get ClinVar field with compatibility for both naming conventions

    Supports both:
    - Old: CLNSIG, CLNDN, CLNSTAR
    - New: ClinVar_CLNSIG, ClinVar_CLNDN, ClinVar_CLNSTAR, ClinVar_CLNREVSTAT
    """
    # Try new naming first (ClinVar_ prefixed)
    new_name = f'ClinVar_{field_name}'
    value = csq_dict.get(new_name, '.')
    if value and value != '.':
        return value

    # Fall back to old naming (without prefix)
    value = csq_dict.get(field_name, '.')
    return value if value and value != '.' else '.'


def evaluate_evoscore_annotation(evo_score):
    """Evaluate EVOScore and return annotation classification

    Classification rules:
    - No value: '.'
    - EVOScore < -11.6: Pathogenic
    - EVOScore > -8.5: Benign
    - Otherwise: Ambiguous
    """
    if not evo_score or evo_score == '.' or evo_score == '':
        return '.'

    try:
        score = float(evo_score)
        if score < -11.6:
            return 'Pathogenic'
        elif score > -8.5:
            return 'Benign'
        else:
            return 'Ambiguous'
    except ValueError:
        return '.'


def evaluate_pangolin_annotation(gain_score, loss_score):
    """Evaluate Pangolin splice scores and return annotation classification

    Classification rules based on the higher score between gain and loss:
    - ≥ 0.5: High likelihood (高度可能导致剪接改变)
    - 0.3 - 0.5: Moderate likelihood (可能影响剪接)
    - 0.1 - 0.3: Low-moderate likelihood (中等可能性)
    - < 0.1: Low likelihood (低可能性)
    - No value: '.'
    """
    # Parse scores
    gain_val = None
    loss_val = None

    # Handle gain score
    if gain_score and gain_score != '.' and gain_score != '':
        try:
            # Handle multiple values (e.g., "0.5&0.3")
            first_gain = gain_score.split('&')[0] if '&' in gain_score else gain_score
            val = float(first_gain)
            # Take absolute value for loss scores (negative indicates splice loss)
            gain_val = abs(val) if val != 0 else None
        except ValueError:
            pass

    # Handle loss score
    if loss_score and loss_score != '.' and loss_score != '':
        try:
            # Handle multiple values (e.g., "-0.88&-0.5")
            first_loss = loss_score.split('&')[0] if '&' in loss_score else loss_score
            val = float(first_loss)
            # Take absolute value for loss scores (negative indicates splice loss)
            loss_val = abs(val) if val != 0 else None
        except ValueError:
            pass

    # Get the maximum score (most significant)
    max_score = None
    if gain_val is not None and loss_val is not None:
        max_score = max(gain_val, loss_val)
    elif gain_val is not None:
        max_score = gain_val
    elif loss_val is not None:
        max_score = loss_val
    else:
        return '.'

    # Classify based on score thresholds
    if max_score >= 0.5:
        return 'High'
    elif max_score >= 0.3:
        return 'Moderate'
    elif max_score >= 0.1:
        return 'Low-mod'
    else:
        return 'Low'


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
    consequence = consequence.replace('_', ' ')
    consequence = consequence.replace('&', ', ')
    return consequence


def get_variant_type(ref, alt):
    """Determine variant type based on REF and ALT"""
    if len(ref) == 1 and len(alt) == 1:
        return 'SNP'
    elif len(ref) == 1 and len(alt) > 1:
        return 'INS'
    elif len(ref) > 1 and len(alt) == 1:
        return 'DEL'
    elif len(ref) > 1 and len(alt) > 1:
        if len(ref) == len(alt):
            return 'MNP'
        else:
            return 'INDEL'
    else:
        return 'OTHER'


def parse_genotype(format_str, sample_str):
    """Parse genotype from FORMAT and sample columns"""
    if not format_str or not sample_str:
        return '.', '.', '.', '.', '.'

    format_fields = format_str.split(':')
    sample_values = sample_str.split(':')

    gt = '.'
    ad = '.'
    dp = '.'
    vaf = '.'
    ps = '.'

    for i, field in enumerate(format_fields):
        if i < len(sample_values):
            if field == 'GT':
                gt = sample_values[i]
            elif field == 'AD':
                ad = sample_values[i]
            elif field == 'DP':
                dp = sample_values[i]
            elif field == 'VAF':
                vaf = sample_values[i]
            elif field == 'PS':
                ps = sample_values[i]

    return gt, ad, dp, vaf, ps


def parse_multiple_samples(format_str, sample_strs):
    """Parse multiple samples from VCF

    Args:
        format_str: FORMAT column string
        sample_strs: List of sample column strings

    Returns:
        List of parsed sample data dictionaries
    """
    samples_data = []
    for sample_str in sample_strs:
        gt, ad, dp, vaf, ps = parse_genotype(format_str, sample_str)
        samples_data.append({
            'gt': gt,
            'ad': ad,
            'dp': dp,
            'vaf': vaf,
            'ps': ps
        })
    return samples_data


def format_hgvs_c(hgvs_c):
    """Remove transcript prefix from HGVS_c

    Example: NM_001005484.2:c.243A>G -> c.243A>G
    """
    if not hgvs_c or hgvs_c == '.':
        return '.'

    # Find the colon and take everything after it
    if ':' in hgvs_c:
        return hgvs_c.split(':', 1)[1]
    return hgvs_c


def format_hgvs_p(hgvs_p):
    """Remove NP prefix from HGVS_p

    Example: NP_001005484.2:p.Ser81= -> p.Ser81=
    """
    if not hgvs_p or hgvs_p == '.':
        return '.'

    # Find the colon and take everything after it
    if ':' in hgvs_p:
        return hgvs_p.split(':', 1)[1]
    return hgvs_p


def determine_zygosity(gt, chrom, sex):
    """Determine zygosity based on genotype, chromosome and sample sex

    Args:
        gt: Genotype string (e.g., '0/0', '0/1', '1/1', '1')
        chrom: Chromosome number/name
        sex: Sample sex ('male' or 'female')

    Returns:
        Zygosity: 'WT', 'Het', 'Hom', or 'Hemi' (for hemizygous)
    """
    if not gt or gt == '.':
        return '.'

    # Normalize genotype
    gt_clean = gt.replace('|', '/')  # Handle phased genotypes

    # Check for male X/Y chromosome hemizygosity
    if sex == 'male':
        # X chromosome in male: hemizygous
        if chrom == 'X':
            if gt_clean in ('0/0', '0', '0|0'):
                return 'WT'
            elif gt_clean in ('1/1', '1', '1|1'):
                return 'Hemi'  # Hemizygous variant
            elif '1' in gt_clean and '0' in gt_clean:
                return 'Hemi'  # Should not normally occur, but handle it
            else:
                return '.'
        # Y chromosome in male: hemizygous
        elif chrom == 'Y':
            if gt_clean in ('0/0', '0', '0|0'):
                return 'WT'
            elif gt_clean in ('1/1', '1', '1|1'):
                return 'Hemi'
            else:
                return '.'

    # Normal diploid chromosomes
    if gt_clean in ('0/0', '0|0', '0'):
        return 'WT'
    elif gt_clean in ('1/1', '1|1', '1'):
        return 'Hom'
    elif '1' in gt_clean and '0' in gt_clean:
        return 'Het'
    else:
        return '.'


def parse_vcf(vcf_path, min_quality=1):
    """Parse VCF file and return list of variant records"""
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

            # Filter out low quality variants (Quality = 0)
            try:
                qual_value = float(qual)
                if qual_value < min_quality:
                    continue
            except ValueError:
                # If quality cannot be parsed as float, skip the record
                if qual == '.' or qual == '0':
                    continue

            format_str = parts[8] if len(parts) > 8 else ''

            # Parse all samples (columns 9 onwards)
            sample_strs = parts[9:] if len(parts) > 9 else []

            info = parse_info_field(info_str)

            # Parse all samples
            samples_data = parse_multiple_samples(format_str, sample_strs)

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
                'csq_entries': csq_entries,
                'samples': samples_data
            }
            records.append(record)

    return records, csq_header


def generate_report(records, output_path, sex='female', sample_names=None):
    """Generate variant report in tab-separated format

    Args:
        records: List of variant records
        output_path: Output file path
        sex: Sample sex ('male' or 'female') for zygosity determination
        sample_names: List of sample names for multi-sample VCF
    """

    # Determine sample names
    if sample_names is None:
        sample_names = ['proband']  # Default single sample name

    # Build headers dynamically based on sample names
    base_headers = [
        'Chromosome',
        'Position',
        'Variant_ID',
        'Ref',
        'Alt',
        'Type',
        'Quality',
        'Filter'
    ]

    # Sample-specific headers
    # First sample: Genotype, Zygosity, PhaseSet, Depth, AD, VAF (no suffix)
    # Other samples: same fields with suffix (e.g., Genotype_father)
    sample_fields = ['Genotype', 'Zygosity', 'PhaseSet', 'Depth', 'AD', 'VAF']

    # First sample headers (no suffix)
    sample_headers = sample_fields.copy()

    # Add headers for additional samples (with suffix)
    for sample_name in sample_names[1:]:
        for field in sample_fields:
            sample_headers.append(f"{field}_{sample_name}")

    # Annotation headers (shared across all samples)
    annotation_headers = [
        'Gene',
        'Transcript',
        'Location',
        'Consequence',
        'Impact',
        'HGVS_c',
        'HGVS_p',
        'Amino_Acids',
        'Cytoband',
        'ClinVar_Sig',
        'ClinVar_RevStat',
        'ClinVar_DN',
        'ClinVar_Star',
        'GnomAD_AF',
        'GnomAD_AF_EAS',
        'GnomAD_nhomalt_XX',
        'GnomAD_nhomalt_XY',
        'Pangolin_Gain',
        'Pangolin_Loss',
        'Pangolin_AN',
        'EVOScore',
        'EVOScore_AN',
        'AlphaMissense_AM',
        'AlphaMissense_AMC',
        'HGNC_ID',
        'Existing_Variation',
        'MAX_AF'
    ]

    headers = base_headers + sample_headers + annotation_headers

    # Sort records by impact, then by chromosome and position
    def sort_key(r):
        # Get first gene symbol from CSQ entries for transcript selection
        first_gene = None
        if r['csq_entries']:
            first_gene = r['csq_entries'][0].get('SYMBOL', None)

        best_csq = get_best_csq_entry(r['csq_entries'], gene=first_gene)
        impact = get_impact_rank(best_csq.get('IMPACT', 'MODIFIER') if best_csq else 'MODIFIER')

        # Sort by impact, then chromosome and position
        try:
            chrom_num = int(r['chrom'])
        except ValueError:
            chrom_num = 100  # X, Y, MT etc
        pos = int(r['pos'])

        return (impact, chrom_num, pos)

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

            # Determine variant type
            var_type = get_variant_type(record['ref'], record['alt'])

            # Build sample data rows
            # First sample: Genotype, Zygosity, PhaseSet, Depth, AD, VAF
            sample_rows = []
            for i, sample_data in enumerate(record['samples']):
                if i < len(sample_names):
                    sample_name = sample_names[i]
                    gt = sample_data['gt']
                    ad = sample_data['ad']
                    dp = sample_data['dp']
                    vaf = sample_data['vaf']
                    ps = sample_data['ps']

                    # Determine zygosity for this sample
                    zygosity = determine_zygosity(gt, record['chrom'], sex)

                    if i == 0:
                        # First sample: no suffix
                        sample_rows.extend([gt, zygosity, ps, dp, ad, vaf])
                    else:
                        # Other samples: with suffix
                        sample_rows.extend([gt, zygosity, ps, dp, ad, vaf])

            # Pad missing samples with '.' if VCF has fewer samples than specified names
            expected_sample_count = len(sample_names)
            actual_sample_count = len(record['samples'])
            if actual_sample_count < expected_sample_count:
                # Each missing sample needs 6 fields
                missing_fields_per_sample = 6
                missing_count = (expected_sample_count - actual_sample_count) * missing_fields_per_sample
                sample_rows.extend(['.'] * missing_count)

            # Gene and transcript info from CSQ
            if best_csq:
                gene = best_csq.get('SYMBOL', '.')
                transcript = best_csq.get('Feature', '.')
                location = extract_transcript_info(best_csq)
                consequence = format_consequence(best_csq.get('Consequence', '.'))
                impact = best_csq.get('IMPACT', '.')
                # Format HGVS to remove transcript/protein prefixes
                hgvs_c = format_hgvs_c(best_csq.get('HGVSc', '.'))
                hgvs_p = format_hgvs_p(best_csq.get('HGVSp', '.'))
                amino_acids = best_csq.get('Amino_acids', '.')
                cytoband = best_csq.get('cytoBand', '.')
                # ClinVar fields with compatibility for both naming conventions
                clinvar_sig = get_clinvar_field(best_csq, 'CLNSIG')
                clinvar_revstat = get_clinvar_field(best_csq, 'CLNREVSTAT')
                clinvar_dn = get_clinvar_field(best_csq, 'CLNDN')
                clinvar_star = get_clinvar_field(best_csq, 'CLNSTAR')
                hgnc_id = best_csq.get('HGNC_ID', '.')
                existing_var = best_csq.get('Existing_variation', '.')
                max_af = best_csq.get('MAX_AF', '.')
                gnomad_af = format_gnomad_af(best_csq)
                gnomad_af_eas = format_gnomad_af_eas(best_csq)
                gnomad_nhomalt_xx = best_csq.get('GnomAD_nhomalt_joint_XX', '.')
                gnomad_nhomalt_xy = best_csq.get('GnomAD_nhomalt_joint_XY', '.')
                pangolin_gain = best_csq.get('Pangolin_gain_score', '.')
                pangolin_loss = best_csq.get('Pangolin_loss_score', '.')
                pangolin_an = evaluate_pangolin_annotation(pangolin_gain, pangolin_loss)
                evo_score = best_csq.get('EVOScore2_EVOScore', '.')
                evo_score_an = evaluate_evoscore_annotation(evo_score)
                alpha_missense_am = best_csq.get('AlphaMissense_AM', '.')
                alpha_missense_amc = best_csq.get('AlphaMissense_AMC', '.')
            else:
                gene = '.'
                transcript = '.'
                location = '.'
                consequence = '.'
                impact = '.'
                hgvs_c = '.'
                hgvs_p = '.'
                amino_acids = '.'
                cytoband = '.'
                clinvar_sig = '.'
                clinvar_revstat = '.'
                clinvar_dn = '.'
                clinvar_star = '.'
                hgnc_id = '.'
                existing_var = '.'
                max_af = '.'
                gnomad_af = '.'
                gnomad_af_eas = '.'
                gnomad_nhomalt_xx = '.'
                gnomad_nhomalt_xy = '.'
                pangolin_gain = '.'
                pangolin_loss = '.'
                pangolin_an = '.'
                evo_score = '.'
                evo_score_an = '.'
                alpha_missense_am = '.'
                alpha_missense_amc = '.'

            row = [
                record['chrom'],
                record['pos'],
                record['var_id'] if record['var_id'] != '.' else '.',
                record['ref'],
                record['alt'],
                var_type,
                record['qual'],
                record['filter']
            ]

            # Add sample data
            row.extend(sample_rows)

            # Add annotation data
            row.extend([
                gene,
                transcript,
                location,
                consequence,
                impact,
                hgvs_c,
                hgvs_p,
                amino_acids,
                cytoband,
                clinvar_sig,
                clinvar_revstat,
                clinvar_dn,
                clinvar_star,
                gnomad_af,
                gnomad_af_eas,
                gnomad_nhomalt_xx,
                gnomad_nhomalt_xy,
                pangolin_gain,
                pangolin_loss,
                pangolin_an,
                evo_score,
                evo_score_an,
                alpha_missense_am,
                alpha_missense_amc,
                hgnc_id,
                existing_var,
                max_af
            ])

            f.write('\t'.join(row) + '\n')

    return len(sorted_records)


def generate_summary_statistics(records, sex='female'):
    """Generate summary statistics for variant report"""
    stats = {
        'total': len(records),
        'pass': 0,
        'filtered': 0,
        'by_type': defaultdict(int),
        'by_impact': defaultdict(int),
        'by_zygosity': defaultdict(int),
        'high_impact_genes': [],
        'moderate_impact_genes': []
    }

    for record in records:
        # Count by filter
        if record['filter'] == 'PASS':
            stats['pass'] += 1
        else:
            stats['filtered'] += 1

        # Count by variant type
        var_type = get_variant_type(record['ref'], record['alt'])
        stats['by_type'][var_type] += 1

        # Count by zygosity (using first sample)
        if record['samples']:
            first_sample_gt = record['samples'][0]['gt']
            zygosity = determine_zygosity(first_sample_gt, record['chrom'], sex)
            stats['by_zygosity'][zygosity] += 1

        # Get first gene symbol for transcript selection
        first_gene = None
        if record['csq_entries']:
            first_gene = record['csq_entries'][0].get('SYMBOL', None)

        # Count by impact
        best_csq = get_best_csq_entry(record['csq_entries'], gene=first_gene)
        if best_csq:
            impact = best_csq.get('IMPACT', 'MODIFIER')
            stats['by_impact'][impact] += 1

            gene = best_csq.get('SYMBOL', '.')
            if impact == 'HIGH' and gene != '.':
                stats['high_impact_genes'].append(gene)
            elif impact == 'MODERATE' and gene != '.':
                stats['moderate_impact_genes'].append(gene)

    return stats


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='VEP VCF Report Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example usage:
  # Single sample
  python vep_report.py -i proband.vep.vcf -o vep_report.txt --sex male

  # Multi-sample trio
  python vep_report.py -i trio.vep.vcf -o trio_report.txt --names proband,father,mother --sex male

  # With custom transcripts file
  python vep_report.py -i proband.vep.vcf -o vep_report.txt -t transcripts.json
        '''
    )

    parser.add_argument('-i', '--input',
                        required=True,
                        help='Input VEP-annotated VCF file')

    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output report file path')

    parser.add_argument('-t', '--transcripts',
                        default=None,
                        help='Transcripts configuration file (default: assets/transcripts.json)')

    parser.add_argument('--sex',
                        choices=['male', 'female'],
                        default='female',
                        help='Sample sex for zygosity determination (default: female). '
                             'Male samples will have X/Y chromosome variants marked as Hemi (hemizygous)')

    parser.add_argument('--names', '-n',
                        default=None,
                        help='Sample names separated by comma (e.g., proband,father,mother). '
                             'First sample\'s data uses base column names (Genotype, Zygosity, etc.), '
                             'other samples use suffixes (Genotype_father, Zygosity_father, etc.)')

    args = parser.parse_args()

    # Parse sample names
    sample_names = None
    if args.names:
        sample_names = [name.strip() for name in args.names.split(',')]
    else:
        sample_names = ['proband']  # Default single sample

    # Set paths
    vcf_path = Path(args.input)
    output_path = Path(args.output)

    # Handle transcripts file
    if args.transcripts:
        global TRANSCRIPTS_FILE
        TRANSCRIPTS_FILE = Path(args.transcripts)

    print(f"Input VCF: {vcf_path}")
    print(f"Output report: {output_path}")
    print(f"Sample names: {', '.join(sample_names)}")
    print(f"Sample sex: {args.sex}")
    if args.transcripts:
        print(f"Transcripts file: {TRANSCRIPTS_FILE}")

    if not vcf_path.exists():
        print(f"Error: VCF file not found: {vcf_path}")
        sys.exit(1)

    # Parse VCF
    print("Parsing VCF file...")
    records, csq_header = parse_vcf(vcf_path)
    print(f"Found {len(records)} variant records (after filtering Quality=0)")

    if not csq_header:
        print("Warning: CSQ header not found, annotation data may be incomplete")

    # Generate report
    print("Generating report...")
    count = generate_report(records, output_path, sex=args.sex, sample_names=sample_names)
    print(f"Report written with {count} entries")

    # Generate summary
    stats = generate_summary_statistics(records, sex=args.sex)
    print("\n=== Variant Summary Statistics ===")
    print(f"Total variants: {stats['total']}")
    print(f"PASS: {stats['pass']}")
    print(f"Filtered: {stats['filtered']}")
    print("\nBy Variant Type:")
    for var_type, cnt in sorted(stats['by_type'].items()):
        print(f"  {var_type}: {cnt}")
    print("\nBy Zygosity:")
    zygosity_order = ['WT', 'Het', 'Hom', 'Hemi', '.']
    for zyg in zygosity_order:
        if zyg in stats['by_zygosity']:
            print(f"  {zyg}: {stats['by_zygosity'][zyg]}")
    for zyg, cnt in sorted(stats['by_zygosity'].items()):
        if zyg not in zygosity_order:
            print(f"  {zyg}: {cnt}")
    print("\nBy Impact:")
    for impact, cnt in sorted(stats['by_impact'].items(), key=lambda x: get_impact_rank(x[0])):
        print(f"  {impact}: {cnt}")
    if stats['high_impact_genes']:
        print(f"\nHIGH impact genes: {', '.join(set(stats['high_impact_genes']))}")
    if stats['moderate_impact_genes']:
        print(f"\nMODERATE impact genes (unique): {len(set(stats['moderate_impact_genes']))}")

    print(f"\nReport saved to: {output_path}")


if __name__ == '__main__':
    main()