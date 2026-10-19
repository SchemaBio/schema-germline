#!/usr/bin/env python3
"""Initialize a single sample VCF using bcftools and a matching FASTA.

Installed in the germline image under /opt/schema-germline/scripts/.
All transformations write into the task work directory, never into input files.
"""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import re
import subprocess
import sys


def run(*args):
    subprocess.run([str(arg) for arg in args], check=True)


def save_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def open_vcf(path):
    with open(path, 'rb') as source:
        compressed = source.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt', encoding='utf-8') if compressed else open(path, encoding='utf-8')


def canonical_chrom(name):
    value = name[3:] if name.startswith('chr') else name
    if value in ('M', 'MT'):
        return 'MT'
    if value in ('X', 'Y') or value in {str(i) for i in range(1, 23)}:
        return value
    return None


def read_fai(reference):
    fai = Path(str(reference) + '.fai')
    if not Path(reference).is_file() or not fai.is_file():
        raise ValueError('Reference FASTA and its existing .fai are required: ' + str(reference))
    contigs = {}
    for line in fai.read_text(encoding='utf-8').splitlines():
        fields = line.split('\t')
        if len(fields) < 5 or not fields[1].isdigit() or int(fields[1]) <= 0:
            raise ValueError('Malformed FASTA index: ' + line)
        if fields[0] in contigs or re.search(r'[\s,<>]', fields[0]):
            raise ValueError('Duplicate or invalid reference contig: ' + fields[0])
        contigs[fields[0]] = int(fields[1])
    if not contigs:
        raise ValueError('Reference FASTA index has no sequences')
    return contigs


def map_contig(name, contigs):
    if name in contigs:
        return name
    canonical = canonical_chrom(name)
    candidates = [contig for contig in contigs if canonical and canonical_chrom(contig) == canonical]
    if len(candidates) != 1:
        raise ValueError(f'Contig {name!r} has no unambiguous match in the reference FASTA')
    return candidates[0]


def parse_gt(fields):
    if len(fields) != 10:
        raise ValueError('Expected exactly one selected sample')
    keys = fields[8].split(':')
    values = fields[9].split(':')
    if 'GT' not in keys or keys.index('GT') >= len(values):
        raise ValueError(f'Missing FORMAT/GT at {fields[0]}:{fields[1]}')
    gt = values[keys.index('GT')]
    parts = re.split(r'[/|]', gt)
    if not parts or any(part in ('', '.') for part in parts):
        return None
    if len(parts) > 2 or any(not part.isdigit() for part in parts):
        raise ValueError(f'Unsupported or malformed GT {gt!r} at {fields[0]}:{fields[1]}')
    alleles = [int(part) for part in parts]
    if any(allele > len(fields[4].split(',')) for allele in alleles):
        raise ValueError(f'GT allele exceeds ALT count at {fields[0]}:{fields[1]}')
    return alleles


def small_variant(ref, alt):
    if not re.fullmatch('[ACGTNacgtn]+', ref) or not re.fullmatch('[ACGTNacgtn]+', alt):
        return False
    return ref.upper() != alt.upper() and (len(ref) != len(alt) or len(ref) == 1)



def validate_config(config):
    if config.get('assembly') not in ('GRCh37', 'GRCh38'):
        raise ValueError('assembly must be GRCh37 or GRCh38; coordinates are not lifted over')
    prefix = config.get('prefix', '')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]{0,127}', prefix):
        raise ValueError('prefix must start with a letter/digit and contain only letters, digits, dot, underscore or hyphen (max 128)')
    if isinstance(config.get('threads', 2), bool) or not isinstance(config.get('threads', 2), int) or config.get('threads', 2) < 1:
        raise ValueError('threads must be a positive integer')


def prepare(config):
    validate_config(config)
    prefix = config['prefix']
    contigs = read_fai(config['reference'])
    # bcftools parses and validates VCF syntax, accepting plain VCF and gzip/BGZF
    # without an input index. Preserve this input copy for line-wise processing.
    run('bcftools', 'view', '--no-version', '-Ov', '-o', 'input.vcf', config['vcf'])
    stats = Counter()
    renames = {}
    declared_lengths = {}
    sample = None
    sample_index = None
    seen_header = False
    with open('input.vcf', encoding='utf-8') as source, open('selected.raw.vcf', 'w', encoding='utf-8') as target:
        for line in source:
            if line.startswith('##contig='):
                name = re.search(r'(?:<|,)ID=([^,>]+)', line)
                length = re.search(r'(?:<|,)length=(\d+)', line)
                if name and length:
                    declared_lengths[name[1]] = int(length[1])
                continue
            if line.startswith('##INFO=<ID=CSQ,'):
                continue
            if line.startswith('##'):
                target.write(line)
                continue
            if line.startswith('#CHROM\t'):
                if seen_header:
                    raise ValueError('Duplicate #CHROM header')
                seen_header = True
                columns = line.rstrip('\r\n').split('\t')
                names = columns[9:]
                if not names or len(set(names)) != len(names):
                    raise ValueError('VCF must contain uniquely named sample columns; sites-only VCF is not supported')
                requested = config.get('sample_name') or None
                if requested is None and len(names) != 1:
                    raise ValueError('Multi-sample VCF requires sample_name; available samples: ' + ', '.join(names))
                sample = requested if requested is not None else names[0]
                if sample not in names:
                    raise ValueError(f'Sample {sample!r} not found; available samples: ' + ', '.join(names))
                if any(char in sample for char in '\r\n\t'):
                    raise ValueError('Invalid sample name')
                sample_index = 9 + names.index(sample)
                for name, length in contigs.items():
                    target.write(f'##contig=<ID={name},length={length}>\n')
                target.write(f'##vcf2anno_assembly={config["assembly"]}\n')
                target.write('\t'.join(columns[:9] + [sample]) + '\n')
                continue
            if not seen_header or not line.strip():
                raise ValueError('VCF record encountered without a valid #CHROM header')
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != len(columns):
                raise ValueError('VCF record column count differs from #CHROM header')
            fields = fields[:9] + [fields[sample_index]]
            stats['input_records'] += 1
            original_chrom = fields[0]
            fields[0] = map_contig(original_chrom, contigs)
            if original_chrom in declared_lengths and declared_lengths[original_chrom] != contigs[fields[0]]:
                raise ValueError(f'Contig length differs from reference: {original_chrom}; check assembly')
            if original_chrom != fields[0]:
                renames[original_chrom] = fields[0]
            gt = parse_gt(fields)
            if gt is None:
                stats['excluded_missing_or_partial_gt'] += 1
                continue
            alt_indices = {allele for allele in gt if allele > 0}
            if not alt_indices:
                stats['excluded_reference_gt'] += 1
                continue
            alts = fields[4].split(',')
            supported = {index for index in alt_indices if small_variant(fields[3], alts[index - 1])}
            stats['unsupported_called_alleles'] += len(alt_indices - supported)
            if not supported:
                stats['excluded_unsupported_type'] += 1
                continue
            # Old VEP annotations must not survive a new annotation run.
            info = [field for field in fields[7].split(';') if field != '.' and not field.startswith('CSQ=')]
            fields[7] = ';'.join(info) or '.'
            target.write('\t'.join(fields) + '\n')
            stats['eligible_input_records'] += 1
    if not seen_header:
        raise ValueError('VCF is missing #CHROM header')
    Path('sample_name.txt').write_text(sample + '\n', encoding='utf-8')
    # Force subset-dependent AC/AN recalculation before splitting; never trust
    # cohort counts to identify the selected sample's alternate alleles.
    run('bcftools', 'view', '--no-version', '--samples-file', 'sample_name.txt',
        '--trim-alt-alleles', '--min-ac', '1:nref', '-Oz', '-o', f'{prefix}.selected.vcf.gz', 'selected.raw.vcf')
    run('bcftools', 'norm', '--no-version', '--check-ref', 'e', '--fasta-ref', config['reference'],
        '--multiallelics', '-any', '--old-rec-tag', 'VCF2ANNO_ORIG', '-Ov', '-o', 'normalized.unsorted.vcf', f'{prefix}.selected.vcf.gz')
    count = 0
    with open('normalized.unsorted.vcf', encoding='utf-8') as source, open('normalized.filtered.vcf', 'w', encoding='utf-8') as target:
        for line in source:
            if line.startswith('#'):
                target.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            gt = parse_gt(fields)
            if gt is None or not any(allele > 0 for allele in gt) or not small_variant(fields[3], fields[4]):
                stats['excluded_after_normalization'] += 1
                continue
            target.write(line)
            count += 1
    normalized = f'{prefix}.normalized.vcf.gz'
    run('bcftools', 'sort', '-Oz', '-o', normalized, 'normalized.filtered.vcf')
    run('bcftools', 'index', '--tbi', '--force', normalized)
    stats['normalized_variants'] = count
    result = {key: stats[key] for key in (
        'input_records', 'eligible_input_records', 'excluded_reference_gt',
        'excluded_missing_or_partial_gt', 'excluded_unsupported_type',
        'unsupported_called_alleles', 'excluded_after_normalization', 'normalized_variants')}
    result.update({'sample_name': sample, 'assembly': config['assembly'], 'chromosome_renames': renames,
                   'quality_filter_applied': False, 'filter_status_filter_applied': False})
    save_json(f'{prefix}.preprocessing.json', result)
    Path('variant_count.txt').write_text(str(count) + '\n', encoding='utf-8')



def main():
    parser = argparse.ArgumentParser(description='Initialize a sample VCF for SNP/InDel annotation')
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--assembly', required=True, choices=['GRCh37', 'GRCh38'])
    parser.add_argument('--reference', required=True)
    parser.add_argument('--sample-name', default=None)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    try:
        prepare(vars(args))
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f'VCF initialization error: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
