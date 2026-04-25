#!/usr/bin/env python3
"""
Mitochondrial VCF Report Generator
Parses VEP-annotated mitochondrial VCF and generates a comprehensive mtDNA variant report

Usage:
    python mt_report.py -i input.vcf -o output.txt [-m mitophen.json]
"""

import re
import sys
import json
import argparse
from pathlib import Path
from collections import defaultdict


# Mitochondrial gene regions (based on rCRS - revised Cambridge Reference Sequence)
# Positions are 1-based inclusive
MT_GENES = {
    # Protein coding genes (13 genes)
    'MT-ND1': {'start': 3307, 'end': 4262, 'type': 'protein', 'desc': 'NADH dehydrogenase subunit 1'},
    'MT-ND2': {'start': 4471, 'end': 5511, 'type': 'protein', 'desc': 'NADH dehydrogenase subunit 2'},
    'MT-CO1': {'start': 5904, 'end': 7445, 'type': 'protein', 'desc': 'Cytochrome c oxidase subunit 1'},
    'MT-CO2': {'start': 7586, 'end': 8269, 'type': 'protein', 'desc': 'Cytochrome c oxidase subunit 2'},
    'MT-ATP8': {'start': 8366, 'end': 8572, 'type': 'protein', 'desc': 'ATP synthase F0 subunit 8'},
    'MT-ATP6': {'start': 8527, 'end': 9207, 'type': 'protein', 'desc': 'ATP synthase F0 subunit 6'},
    'MT-CO3': {'start': 9207, 'end': 9990, 'type': 'protein', 'desc': 'Cytochrome c oxidase subunit 3'},
    'MT-ND3': {'start': 10059, 'end': 10404, 'type': 'protein', 'desc': 'NADH dehydrogenase subunit 3'},
    'MT-ND4L': {'start': 10470, 'end': 10766, 'type': 'protein', 'desc': 'NADH dehydrogenase subunit 4L'},
    'MT-ND4': {'start': 10760, 'end': 12137, 'type': 'protein', 'desc': 'NADH dehydrogenase subunit 4'},
    'MT-ND5': {'start': 12337, 'end': 14148, 'type': 'protein', 'desc': 'NADH dehydrogenase subunit 5'},
    'MT-ND6': {'start': 14149, 'end': 14673, 'type': 'protein', 'desc': 'NADH dehydrogenase subunit 6', 'strand': -1},
    'MT-CYB': {'start': 14747, 'end': 15887, 'type': 'protein', 'desc': 'Cytochrome b'},
    # rRNA genes (2 genes)
    'MT-RNR1': {'start': 648, 'end': 1601, 'type': 'rRNA', 'desc': '12S ribosomal RNA'},
    'MT-RNR2': {'start': 1671, 'end': 3229, 'type': 'rRNA', 'desc': '16S ribosomal RNA'},
    # tRNA genes (22 genes)
    'MT-TF': {'start': 577, 'end': 647, 'type': 'tRNA', 'desc': 'tRNA-Phe'},
    'MT-TV': {'start': 1602, 'end': 1670, 'type': 'tRNA', 'desc': 'tRNA-Val'},
    'MT-TL1': {'start': 3230, 'end': 3304, 'type': 'tRNA', 'desc': 'tRNA-Leu(UUR)'},
    'MT-TI': {'start': 4263, 'end': 4331, 'type': 'tRNA', 'desc': 'tRNA-Ile'},
    'MT-TQ': {'start': 4332, 'end': 4400, 'type': 'tRNA', 'desc': 'tRNA-Met'},
    'MT-TM': {'start': 4402, 'end': 4469, 'type': 'tRNA', 'desc': 'tRNA-Met'},
    'MT-TW': {'start': 5512, 'end': 5579, 'type': 'tRNA', 'desc': 'tRNA-Trp'},
    'MT-TA': {'start': 5587, 'end': 5655, 'type': 'tRNA', 'desc': 'tRNA-Ala'},
    'MT-TN': {'start': 5657, 'end': 5729, 'type': 'tRNA', 'desc': 'tRNA-Asn'},
    'MT-TC': {'start': 5826, 'end': 5891, 'type': 'tRNA', 'desc': 'tRNA-Cys'},
    'MT-TY': {'start': 5892, 'end': 5900, 'type': 'tRNA', 'desc': 'tRNA-Tyr'},
    'MT-TS1': {'start': 7446, 'end': 7514, 'type': 'tRNA', 'desc': 'tRNA-Ser(UCN)'},
    'MT-TD': {'start': 7518, 'end': 7585, 'type': 'tRNA', 'desc': 'tRNA-Asp'},
    'MT-TK': {'start': 8295, 'end': 8364, 'type': 'tRNA', 'desc': 'tRNA-Lys'},
    'MT-TG': {'start': 9991, 'end': 10058, 'type': 'tRNA', 'desc': 'tRNA-Gly'},
    'MT-TR': {'start': 10405, 'end': 10469, 'type': 'tRNA', 'desc': 'tRNA-Arg'},
    'MT-TH': {'start': 12138, 'end': 12206, 'type': 'tRNA', 'desc': 'tRNA-His'},
    'MT-TS2': {'start': 12207, 'end': 12272, 'type': 'tRNA', 'desc': 'tRNA-Ser(AGY)'},
    'MT-TL2': {'start': 12273, 'end': 12336, 'type': 'tRNA', 'desc': 'tRNA-Leu(CUN)'},
    'MT-TE': {'start': 14674, 'end': 14742, 'type': 'tRNA', 'desc': 'tRNA-Glu'},
    'MT-TT': {'start': 15888, 'end': 15953, 'type': 'tRNA', 'desc': 'tRNA-Thr'},
    'MT-TP': {'start': 15955, 'end': 16023, 'type': 'tRNA', 'desc': 'tRNA-Pro'},
    # D-loop / Control region (non-coding)
    'DLOOP': {'start': 16024, 'end': 16569, 'type': 'D-loop', 'desc': 'Control region / D-loop'},
}

# Non-coding regions (intergenic)
MT_NONCODING = {
    'NC1': {'start': 1, 'end': 576, 'desc': 'Non-coding region before MT-TF'},
    'NC2': {'start': 5756, 'end': 5825, 'desc': 'Non-coding region between MT-TN and MT-TC'},
}


def get_mt_gene_by_position(pos):
    """Determine which mitochondrial gene a position falls into

    Args:
        pos: Position on MT chromosome (1-based)

    Returns:
        Gene name and gene type, or 'Intergenic' if not in a gene
    """
    try:
        pos_int = int(pos)
    except ValueError:
        return 'Intergenic', 'unknown'

    # Check all gene regions
    for gene_name, gene_info in MT_GENES.items():
        if gene_info['start'] <= pos_int <= gene_info['end']:
            return gene_name, gene_info['type']

    # Check non-coding regions
    for region_name, region_info in MT_NONCODING.items():
        if region_info['start'] <= pos_int <= region_info['end']:
            return 'Intergenic', 'noncoding'

    return 'Intergenic', 'unknown'


# Mitophen database cache
_mitophen_cache = None


def load_mitophen(mitophen_path):
    """Load mitophen database for phenotype annotation"""
    global _mitophen_cache
    if _mitophen_cache is None and mitophen_path:
        try:
            with open(mitophen_path, 'r', encoding='utf-8') as f:
                _mitophen_cache = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            print(f"Warning: Could not load mitophen database: {e}")
            _mitophen_cache = {}
    return _mitophen_cache


def get_mitophen_phenotypes(pos, ref, alt, mitophen_data):
    """Get phenotypes for a mitochondrial variant from mitophen database

    Args:
        pos: Position on MT chromosome
        ref: Reference allele
        alt: Alternate allele
        mitophen_data: Loaded mitophen database

    Returns:
        Tuple of (variant_name, phenotypes_str) or ('.', '.') if not found
    """
    if not mitophen_data:
        return '.', '.'

    pos_map = mitophen_data.get('position_map', {})
    # Convert pos to string for lookup (position_map keys are strings)
    pos_str = str(int(pos))

    # Check if position has known variants
    if pos_str not in pos_map:
        return '.', '.'

    # Try to match variant by ref/alt
    for variant_info in pos_map[pos_str]:
        variant_name = variant_info.get('name', '')
        # Parse variant name like "m.583G>A" to extract ref and alt
        match = re.match(r'm\.(\d+)([ACGT])>([ACGT])', variant_name)
        if match:
            v_ref = match.group(2)
            v_alt = match.group(3)
            if v_ref == ref and v_alt == alt:
                phenotypes = variant_info.get('phenotypes', [])
                # Get top phenotypes (skip "All" and generic terms)
                significant_phenotypes = []
                for ph in phenotypes:
                    if ph['code'] not in ['HP:0000001', 'HP:0000118', 'HP:0000005']:
                        significant_phenotypes.append(ph['name'])
                # Return top 10 phenotypes
                phenotype_str = '; '.join(significant_phenotypes[:10]) if significant_phenotypes else '.'
                return variant_name, phenotype_str

    return '.', '.'


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


def parse_genotype(format_str, sample_str):
    """Parse genotype from FORMAT and sample columns"""
    if not format_str or not sample_str:
        return '.', '.', '.', '.'

    format_fields = format_str.split(':')
    sample_values = sample_str.split(':')

    gt = '.'
    ad = '.'
    af = '.'
    dp = '.'

    for i, field in enumerate(format_fields):
        if i < len(sample_values):
            if field == 'GT':
                gt = sample_values[i]
            elif field == 'AD':
                ad = sample_values[i]
            elif field == 'AF':
                af = sample_values[i]
            elif field == 'DP':
                dp = sample_values[i]

    return gt, ad, af, dp


def get_impact_rank(impact):
    """Rank impact for sorting (HIGH=0, MODERATE=1, LOW=2, MODIFIER=3)"""
    impact_order = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
    return impact_order.get(impact, 4)


def format_hgvs_c(hgvs_c):
    """Remove transcript prefix from HGVS_c"""
    if not hgvs_c or hgvs_c == '.':
        return '.'
    if ':' in hgvs_c:
        return hgvs_c.split(':', 1)[1]
    return hgvs_c


def format_hgvs_p(hgvs_p):
    """Remove protein prefix from HGVS_p"""
    if not hgvs_p or hgvs_p == '.':
        return '.'
    if ':' in hgvs_p:
        return hgvs_p.split(':', 1)[1]
    return hgvs_p


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


def get_best_csq_entry(csq_entries):
    """Get the most significant CSQ entry (highest impact)"""
    if not csq_entries:
        return None

    # Sort by impact rank
    sorted_entries = sorted(csq_entries, key=lambda x: get_impact_rank(x.get('IMPACT', 'MODIFIER')))
    return sorted_entries[0]


def format_consequence(consequence):
    """Format consequence for readability"""
    consequence = consequence.replace('_', ' ')
    consequence = consequence.replace('&', ', ')
    return consequence


def format_mt_hgvs(hgvsg):
    """Format mitochondrial HGVS notation

    Example: MT:m.23T>G -> m.23T>G
    """
    if not hgvsg or hgvsg == '.':
        return '.'
    if ':' in hgvsg:
        return hgvsg.split(':', 1)[1]
    return hgvsg


def extract_rsid(existing_var):
    """Extract rsid from Existing_variation field

    Example: 'rs2124593767&COSV62293932' -> 'rs2124593767'
    Example: 'rs1556422421' -> 'rs1556422421'
    Example: '.' -> '.'
    """
    if not existing_var or existing_var == '.':
        return '.'

    # Find all rsids using regex
    rsids = re.findall(r'rs\d+', existing_var)

    if rsids:
        # Return first rsid (or join multiple with comma)
        return rsids[0] if len(rsids) == 1 else ','.join(rsids)

    return '.'


def calculate_heteroplasmy(ad, af):
    """Calculate heteroplasmy level from AD or AF

    For mitochondrial DNA, heteroplasmy is the proportion of mutant alleles
    """
    # Prefer AF field directly if available
    if af and af != '.' and af != '':
        try:
            return float(af)
        except ValueError:
            pass

    # Calculate from AD if AF not available
    if ad and ad != '.' and ad != '':
        try:
            ad_values = ad.split(',')
            if len(ad_values) >= 2:
                ref_depth = int(ad_values[0])
                alt_depths = [int(x) for x in ad_values[1:] if x]
                total_alt = sum(alt_depths)
                total = ref_depth + total_alt
                if total > 0:
                    return total_alt / total
        except ValueError:
            pass

    return 0.0


def classify_heteroplasmy(heteroplasmy):
    """Classify heteroplasmy level

    Classification:
    - >= 0.95: Homoplasmy
    - 0.60 - 0.95: High heteroplasmy
    - 0.20 - 0.60: Moderate heteroplasmy
    - 0.01 - 0.20: Low heteroplasmy
    - < 0.01: Very low heteroplasmy
    """
    if heteroplasmy >= 0.95:
        return 'Homoplasmy'
    elif heteroplasmy >= 0.60:
        return 'High'
    elif heteroplasmy >= 0.20:
        return 'Moderate'
    elif heteroplasmy >= 0.01:
        return 'Low'
    else:
        return 'Very-low'


def parse_vcf(vcf_path):
    """Parse VCF file and return list of mitochondrial variant records"""
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

            # Only process mitochondrial variants
            if chrom != 'MT':
                continue

            format_str = parts[8] if len(parts) > 8 else ''
            sample_str = parts[9] if len(parts) > 9 else ''

            info = parse_info_field(info_str)

            # Parse genotype
            gt, ad, af, dp = parse_genotype(format_str, sample_str)

            # Calculate heteroplasmy
            heteroplasmy = calculate_heteroplasmy(ad, af)
            heteroplasmy_class = classify_heteroplasmy(heteroplasmy)

            # Parse CSQ
            csq_str = info.get('CSQ', '')
            csq_entries = parse_csq_field(csq_str, csq_header) if csq_header else []

            # Get additional INFO fields
            tlod = info.get('TLOD', '.')
            popaf = info.get('POPAF', '.')
            germq = info.get('GERMQ', '.')
            strandq = info.get('STRANDQ', '.')
            contq = info.get('CONTQ', '.')
            seqq = info.get('SEQQ', '.')
            mbq = info.get('MBQ', '.')
            mmq = info.get('MMQ', '.')
            mfrl = info.get('MFRL', '.')

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
                'gt': gt,
                'ad': ad,
                'af': af,
                'dp': dp,
                'heteroplasmy': heteroplasmy,
                'heteroplasmy_class': heteroplasmy_class,
                'tlod': tlod,
                'popaf': popaf,
                'germq': germq,
                'strandq': strandq,
                'contq': contq,
                'seqq': seqq,
                'mbq': mbq,
                'mmq': mmq,
                'mfrl': mfrl
            }
            records.append(record)

    return records, csq_header


def generate_report(records, output_path, mitophen_data=None):
    """Generate mitochondrial variant report in tab-separated format"""

    # Define column headers
    headers = [
        'Chromosome',
        'Position',
        'MT_Gene',
        'MT_Gene_Type',
        'Mitophen_Variant',
        'Mitophen_Phenotypes',
        'MT_HGVS',
        'Ref',
        'Alt',
        'Type',
        'Filter',
        'Genotype',
        'Heteroplasmy',
        'Heteroplasmy_Class',
        'Depth',
        'AD',
        'AF',
        'Gene',
        'Feature',
        'Consequence',
        'Impact',
        'HGVS_c',
        'HGVS_p',
        'Amino_Acids',
        'Protein_Position',
        'ClinVar_Sig',
        'ClinVar_DN',
        'ClinVar_Star',
        'GnomAD_AF',
        'GnomAD_AF_EAS',
        'dbSNP',
        'MAX_AF',
        'TLOD',
        'POPAF',
        'GERMQ',
        'STRANDQ',
        'CONTQ',
        'SEQQ',
        'MBQ',
        'MMQ',
        'MFRL'
    ]

    # Sort records by heteroplasmy (descending), then by position
    sorted_records = sorted(records, key=lambda r: (-r['heteroplasmy'], int(r['pos'])))

    # Write report
    with open(output_path, 'w', encoding='utf-8') as f:
        # Write header
        f.write('\t'.join(headers) + '\n')

        # Write data rows
        for record in sorted_records:
            info = record['info']

            best_csq = get_best_csq_entry(record['csq_entries'])

            # Determine MT gene by position
            mt_gene, mt_gene_type = get_mt_gene_by_position(record['pos'])

            # Get Mitophen phenotypes
            mitophen_var, mitophen_pheno = get_mitophen_phenotypes(
                record['pos'], record['ref'], record['alt'], mitophen_data
            )

            # Determine variant type
            var_type = get_variant_type(record['ref'], record['alt'])

            # Get MT HGVS notation
            mt_hgvs = '.'
            if best_csq:
                mt_hgvs = format_mt_hgvs(best_csq.get('HGVSg', '.'))

            # Format heteroplasmy
            heteroplasmy_str = f"{record['heteroplasmy']:.4f}" if record['heteroplasmy'] > 0 else '0'

            # Gene and transcript info from CSQ
            if best_csq:
                gene = best_csq.get('SYMBOL', '.')
                feature = best_csq.get('Feature', '.')
                consequence = format_consequence(best_csq.get('Consequence', '.'))
                impact = best_csq.get('IMPACT', '.')
                hgvs_c = format_hgvs_c(best_csq.get('HGVSc', '.'))
                hgvs_p = format_hgvs_p(best_csq.get('HGVSp', '.'))
                amino_acids = best_csq.get('Amino_acids', '.')
                protein_position = best_csq.get('Protein_position', '.')
                # ClinVar fields
                clinvar_sig = best_csq.get('CLNSIG', '.') if best_csq.get('CLNSIG', '.') else '.'
                clinvar_dn = best_csq.get('CLNDN', '.') if best_csq.get('CLNDN', '.') else '.'
                clinvar_star = best_csq.get('CLNSTAR', '.') if best_csq.get('CLNSTAR', '.') else '.'
                # GnomAD
                gnomad_af = best_csq.get('GnomAD_AF_joint', '.') if best_csq.get('GnomAD_AF_joint', '.') else '.'
                gnomad_af_eas = best_csq.get('GnomAD_AF_joint_eas', '.') if best_csq.get('GnomAD_AF_joint_eas', '.') else '.'
                existing_var_raw = best_csq.get('Existing_variation', '.')
                existing_var = extract_rsid(existing_var_raw)
                max_af = best_csq.get('MAX_AF', '.')
            else:
                gene = '.'
                feature = '.'
                consequence = 'intergenic_variant'
                impact = 'MODIFIER'
                hgvs_c = '.'
                hgvs_p = '.'
                amino_acids = '.'
                protein_position = '.'
                clinvar_sig = '.'
                clinvar_dn = '.'
                clinvar_star = '.'
                gnomad_af = '.'
                gnomad_af_eas = '.'
                existing_var = '.'
                max_af = '.'

            row = [
                record['chrom'],
                record['pos'],
                mt_gene,
                mt_gene_type,
                mitophen_var,
                mitophen_pheno,
                mt_hgvs,
                record['ref'],
                record['alt'],
                var_type,
                record['filter'],
                record['gt'],
                heteroplasmy_str,
                record['heteroplasmy_class'],
                record['dp'],
                record['ad'],
                record['af'],
                gene,
                feature,
                consequence,
                impact,
                hgvs_c,
                hgvs_p,
                amino_acids,
                protein_position,
                clinvar_sig,
                clinvar_dn,
                clinvar_star,
                gnomad_af,
                gnomad_af_eas,
                existing_var,
                max_af,
                record['tlod'],
                record['popaf'],
                record['germq'],
                record['strandq'],
                record['contq'],
                record['seqq'],
                record['mbq'],
                record['mmq'],
                record['mfrl']
            ]

            f.write('\t'.join(row) + '\n')

    return len(sorted_records)


def generate_summary_statistics(records):
    """Generate summary statistics for mitochondrial variant report"""
    stats = {
        'total': len(records),
        'by_type': defaultdict(int),
        'by_impact': defaultdict(int),
        'by_heteroplasmy': defaultdict(int),
        'by_gene_type': defaultdict(int),
        'homoplasmy_count': 0,
        'high_heteroplasmy_count': 0,
        'gene_variants': defaultdict(int)
    }

    for record in records:
        # Count by variant type
        var_type = get_variant_type(record['ref'], record['alt'])
        stats['by_type'][var_type] += 1

        # Count by heteroplasmy class
        stats['by_heteroplasmy'][record['heteroplasmy_class']] += 1

        if record['heteroplasmy'] >= 0.95:
            stats['homoplasmy_count'] += 1
        elif record['heteroplasmy'] >= 0.60:
            stats['high_heteroplasmy_count'] += 1

        # Count by MT gene type
        mt_gene, mt_gene_type = get_mt_gene_by_position(record['pos'])
        stats['by_gene_type'][mt_gene_type] += 1

        # Count by impact
        best_csq = get_best_csq_entry(record['csq_entries'])
        if best_csq:
            impact = best_csq.get('IMPACT', 'MODIFIER')
            stats['by_impact'][impact] += 1

            gene = best_csq.get('SYMBOL', '.')
            if gene and gene != '.':
                stats['gene_variants'][gene] += 1

    return stats


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Mitochondrial VCF Report Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example usage:
  python mt_report.py -i proband.mt.vep.vcf -o mt_report.txt
        '''
    )

    parser.add_argument('-i', '--input',
                        required=True,
                        help='Input VEP-annotated mitochondrial VCF file')

    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output report file path')

    parser.add_argument('-m', '--mitophen',
                        default=None,
                        help='Mitophen database JSON file for phenotype annotation')

    args = parser.parse_args()

    # Set paths
    vcf_path = Path(args.input)
    output_path = Path(args.output)
    mitophen_path = args.mitophen

    print(f"Input VCF: {vcf_path}")
    print(f"Output report: {output_path}")
    if mitophen_path:
        print(f"Mitophen database: {mitophen_path}")

    if not vcf_path.exists():
        print(f"Error: VCF file not found: {vcf_path}")
        sys.exit(1)

    # Load mitophen database
    mitophen_data = load_mitophen(mitophen_path)
    if mitophen_data:
        metadata = mitophen_data.get('metadata', {})
        print(f"Loaded mitophen v{metadata.get('version', '?')} with {metadata.get('variants_with_phenotypes', 0)} annotated variants")

    # Parse VCF
    print("Parsing VCF file...")
    records, csq_header = parse_vcf(vcf_path)
    print(f"Found {len(records)} mitochondrial variant records")

    if not csq_header:
        print("Warning: CSQ header not found, annotation data may be incomplete")

    # Generate report
    print("Generating report...")
    count = generate_report(records, output_path, mitophen_data)
    print(f"Report written with {count} entries")

    # Generate summary
    stats = generate_summary_statistics(records)
    print("\n=== Mitochondrial Variant Summary Statistics ===")
    print(f"Total variants: {stats['total']}")
    print(f"Homoplasmy (≥95%): {stats['homoplasmy_count']}")
    print(f"High heteroplasmy (60-95%): {stats['high_heteroplasmy_count']}")

    print("\nBy Variant Type:")
    for var_type, cnt in sorted(stats['by_type'].items()):
        print(f"  {var_type}: {cnt}")

    print("\nBy Heteroplasmy Level:")
    het_order = ['Homoplasmy', 'High', 'Moderate', 'Low', 'Very-low']
    for het in het_order:
        if het in stats['by_heteroplasmy']:
            print(f"  {het}: {stats['by_heteroplasmy'][het]}")

    print("\nBy Impact:")
    for impact, cnt in sorted(stats['by_impact'].items(), key=lambda x: get_impact_rank(x[0])):
        print(f"  {impact}: {cnt}")

    print("\nBy MT Gene Type:")
    gene_type_order = ['protein', 'rRNA', 'tRNA', 'D-loop', 'noncoding', 'unknown']
    for gt in gene_type_order:
        if gt in stats['by_gene_type']:
            print(f"  {gt}: {stats['by_gene_type'][gt]}")

    if stats['gene_variants']:
        print("\nVariants by Gene:")
        for gene, cnt in sorted(stats['gene_variants'].items()):
            print(f"  {gene}: {cnt}")

    print(f"\nReport saved to: {output_path}")


if __name__ == '__main__':
    main()