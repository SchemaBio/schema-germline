#!/usr/bin/env python3
"""Infer a sample's sex from X variants before multiallelic splitting.

Both scripts are installed together in the germline image. No script is a WDL
File input; the workflow invokes the fixed /opt/schema-germline/scripts/ path.
"""
import argparse
import json
import math
from pathlib import Path
import sys

from initialize_vcf import canonical_chrom, open_vcf, parse_gt, small_variant, save_json

SEX_WARNING = "由 VCF 的 X 染色体位点比例推算，可能受 X 染色体变异、检测区域及变异检出方式影响。"


def infer_sex(x_sum, x_hom, x_haploid, minimum=100, female_max=0.5, male_min=0.7):
    if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 1:
        raise ValueError('sex_min_sites must be a positive integer')
    if not (math.isfinite(female_max) and math.isfinite(male_min) and 0 <= female_max < male_min <= 1):
        raise ValueError('Require 0 <= female_max_ratio < male_min_ratio <= 1')
    ratio = (x_hom + x_haploid) / x_sum if x_sum else None
    if x_sum == 0:
        sex, reason = 'unknown', 'no_x_variants'
    elif x_sum < minimum:
        sex, reason = 'unknown', 'insufficient_x_variants'
    elif ratio <= female_max:
        sex, reason = 'female', 'ratio_at_or_below_female_max'
    elif ratio >= male_min:
        sex, reason = 'male', 'ratio_at_or_above_male_min'
    else:
        sex, reason = 'unknown', 'ratio_in_uncertain_interval'
    return {
        'method': 'x_variant_homozygosity', 'sex': sex, 'reason': reason,
        'x_sum': x_sum, 'x_hom': x_hom, 'x_haploid_alt': x_haploid,
        'numerator': x_hom + x_haploid, 'ratio': ratio,
        'minimum_sites': minimum, 'female_max_ratio': female_max, 'male_min_ratio': male_min,
        'counting_basis': 'selected-sample original variant records before multiallelic splitting; PAR included',
        'warning': SEX_WARNING,
    }



def infer_from_vcf(vcf, minimum=100, female_max=0.5, male_min=0.7):
    # Validate even empty VCFs before inspecting any records.
    infer_sex(0, 0, 0, minimum, female_max, male_min)
    x_sum = x_hom = x_haploid = 0
    with open_vcf(vcf) as source:
        for line in source:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if canonical_chrom(fields[0]) != 'X':
                continue
            gt = parse_gt(fields)
            if gt is None:
                continue
            alts = fields[4].split(',')
            supported = {allele for allele in gt if allele > 0 and small_variant(fields[3], alts[allele - 1])}
            if not supported:
                continue
            x_sum += 1
            if len(gt) == 2 and gt[0] == gt[1] and gt[0] in supported:
                x_hom += 1
            elif len(gt) == 1 and gt[0] in supported:
                x_haploid += 1
    return infer_sex(x_sum, x_hom, x_haploid, minimum, female_max, male_min)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', required=True, help='Selected sample VCF before multiallelic splitting')
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--sample-name', required=True)
    parser.add_argument('--minimum-sites', type=int, default=100)
    parser.add_argument('--female-max-ratio', type=float, default=0.5)
    parser.add_argument('--male-min-ratio', type=float, default=0.7)
    args = parser.parse_args()
    try:
        result = infer_from_vcf(args.vcf, args.minimum_sites, args.female_max_ratio, args.male_min_ratio)
        result['sample_name'] = args.sample_name
        save_json(f'{args.prefix}.sex_inference.json', result)
        ratio = 'NA' if result['ratio'] is None else f'{result["ratio"]:.6f}'
        lines = [f'Sample: {args.sample_name}', f'Inferred sex: {result["sex"]}', f'Reason: {result["reason"]}',
                 f'X(sum): {result["x_sum"]}', f'X(hom): {result["x_hom"]}', f'X(haploid ALT): {result["x_haploid_alt"]}',
                 f'Ratio (X(hom) + X(haploid ALT)) / X(sum): {ratio}',
                 f'Minimum sites: {result["minimum_sites"]}',
                 f'Female <= {result["female_max_ratio"]}; male >= {result["male_min_ratio"]}', SEX_WARNING]
        Path(f'{args.prefix}.sex_inference.txt').write_text('\n'.join(lines)+'\n', encoding='utf-8')
        Path('inferred_sex.txt').write_text(result['sex']+'\n', encoding='utf-8')
    except (ValueError, OSError) as error:
        print(f'VCF sex inference error: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
