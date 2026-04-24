#!/usr/bin/env python3
"""
全外显子质控报表生成脚本
从test文件夹的中间文件提取关键质控信息，生成JSON格式报表
"""

import json
import argparse


def parse_fastp_stats(fastp_file: str) -> dict:
    """
    从fastp_stats.json提取过滤前后的关键指标

    返回:
        dict: 包含before_filtering和after_filtering的质控数据
    """
    with open(fastp_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    summary = data.get('summary', {})

    result = {}

    # 提取过滤前数据
    before = summary.get('before_filtering', {})
    result['before_filtering'] = {
        'total_reads': before.get('total_reads', 0),
        'total_bases': before.get('total_bases', 0),
        'q20_rate': before.get('q20_rate', 0),
        'q30_rate': before.get('q30_rate', 0),
        'read1_mean_length': before.get('read1_mean_length', 0),
        'read2_mean_length': before.get('read2_mean_length', 0),
        'gc_content': before.get('gc_content', 0)
    }

    # 提取过滤后数据
    after = summary.get('after_filtering', {})
    result['after_filtering'] = {
        'total_reads': after.get('total_reads', 0),
        'total_bases': after.get('total_bases', 0),
        'q20_rate': after.get('q20_rate', 0),
        'q30_rate': after.get('q30_rate', 0),
        'read1_mean_length': after.get('read1_mean_length', 0),
        'read2_mean_length': after.get('read2_mean_length', 0),
        'gc_content': after.get('gc_content', 0)
    }

    return result


def parse_xamdst_report(xamdst_file: str) -> dict:
    """
    从xamdst.report.json提取比对和覆盖度关键指标

    返回:
        dict: 包含比对统计和覆盖度数据
    """
    with open(xamdst_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    result = {}

    # 提取比对统计
    total = data.get('total', {})
    result['mapped_reads'] = total.get('mapped_reads', 0)
    result['mapped_reads_fraction'] = total.get('mapped_reads_fraction', 0)

    # 提取插入片段大小
    insert_size = data.get('insert_size', {})
    result['insert_size_average'] = insert_size.get('average', 0)
    result['insert_size_median'] = insert_size.get('median', 0)

    # 提取目标区域统计
    target = data.get('target', {})
    result['target_data_fraction_all'] = target.get('target_data_fraction_all', 0)
    result['average_depth'] = target.get('average_depth', 0)
    result['average_depth_rmdup'] = target.get('average_depth_rmdup', 0)
    result['region_length'] = target.get('region_length', 0)

    # 提取覆盖度指标（去重后数据）
    coverage_rmdup = target.get('coverage_rmdup', {})
    result['coverage_gt_0x'] = coverage_rmdup.get('gt_0x', 0)
    result['coverage_gte_30x'] = coverage_rmdup.get('gte_30x', 0)
    result['coverage_gte_100x'] = coverage_rmdup.get('gte_100x', 0)

    # gt_0_2_avg 使用去重前数据
    coverage = target.get('coverage', {})
    result['coverage_gt_0_2_avg'] = coverage.get('gt_0_2_avg', 0)

    return result


def parse_mt_xamdst_report(mt_xamdst_file: str) -> dict:
    """
    从mt.xamdst.report.json提取线粒体覆盖度指标

    返回:
        dict: 包含线粒体去重后覆盖度数据
    """
    with open(mt_xamdst_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    result = {}

    target = data.get('target', {})

    result['mt_average_depth'] = target.get('average_depth', 0)
    result['mt_average_depth_rmdup'] = target.get('average_depth_rmdup', 0)

    coverage_rmdup = target.get('coverage_rmdup', {})

    result['mt_coverage_gt_0x'] = coverage_rmdup.get('gt_0x', 0)

    # gte_1000x 在 custom 字段中
    custom = coverage_rmdup.get('custom', {})
    result['mt_coverage_gte_1000x'] = custom.get('gte_1000x', 0)

    return result


def parse_fingerprint(fingerprint_file: str) -> dict:
    """
    从fingerprint.json提取指纹信息

    返回:
        dict: 包含fingerprint和fingerprint_hash
    """
    with open(fingerprint_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    result = {
        'fingerprint': data.get('fingerprint', ''),
        'fingerprint_hash': data.get('fingerprint_hash', '')
    }

    return result


def parse_metrics(metrics_file: str) -> dict:
    """
    从metrics.txt提取Picard比对质量指标

    返回:
        dict: 包含PF_MISMATCH_RATE
    """
    result = {}

    with open(metrics_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # 找到METRICS CLASS行和数据行
    for i, line in enumerate(lines):
        if line.startswith('## METRICS CLASS'):
            # 下一行是列名，再下一行开始是数据
            header_line = lines[i + 1].strip()
            data_lines = []
            for j in range(i + 2, len(lines)):
                data_line = lines[j].strip()
                if data_line and not data_line.startswith('##'):
                    data_lines.append(data_line)
                else:
                    break

            # 解析列名
            headers = header_line.split('\t')

            # 找PAIR行的数据（整体统计）
            for data_line in data_lines:
                values = data_line.split('\t')
                if values[0] == 'PAIR':
                    # 找到PF_MISMATCH_RATE列
                    for idx, header in enumerate(headers):
                        if header == 'PF_MISMATCH_RATE':
                            result['pf_mismatch_rate'] = float(values[idx]) if idx < len(values) else 0
                            break
                    break
            break

    return result


def parse_hs_metrics(hs_file: str) -> dict:
    """
    从hs.txt提取Picard HsMetrics靶向测序质控指标

    返回:
        dict: 包含覆盖度和富集效率相关指标
    """
    result = {}

    with open(hs_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # 找到METRICS CLASS行和数据行
    for i, line in enumerate(lines):
        if line.startswith('## METRICS CLASS') and 'HsMetrics' in line:
            # 下一行是列名，再下一行开始是数据
            header_line = lines[i + 1].strip()
            data_line = lines[i + 2].strip()

            if data_line and not data_line.startswith('##'):
                headers = header_line.split('\t')
                values = data_line.split('\t')

                # 提取关键指标
                key_fields = [
                    'MEAN_TARGET_COVERAGE',
                    'MEDIAN_TARGET_COVERAGE',
                    'PCT_TARGET_BASES_30X',
                    'PCT_TARGET_BASES_100X',
                    'ZERO_CVG_TARGETS_PCT',
                    'FOLD_ENRICHMENT',
                    'HET_SNP_SENSITIVITY',
                    'FOLD_80_BASE_PENALTY'
                ]

                for idx, header in enumerate(headers):
                    if header in key_fields and idx < len(values):
                        key_name = header.lower()
                        try:
                            result[key_name] = float(values[idx])
                        except ValueError:
                            result[key_name] = values[idx]
            break

    return result


def parse_sry_count(sry_file: str, cutoff: int) -> dict:
    """
    从SRY.count.txt提取SRY计数并判断性别

    参数:
        cutoff: SRY计数阈值，大于此值为male

    返回:
        dict: 包含sry_count和predicted_gender
    """
    result = {}

    with open(sry_file, 'r', encoding='utf-8') as f:
        content = f.read().strip()

    try:
        sry_count = int(content)
    except ValueError:
        sry_count = 0

    result['sry_count'] = sry_count
    result['predicted_gender'] = 'male' if sry_count > cutoff else 'female'

    return result


def parse_library_complexity(lc_file: str) -> dict:
    """
    从GATK EstimateLibraryComplexity输出文件提取文库复杂度指标

    返回:
        dict: 包含estimated_library_size, percent_duplication等
    """
    result = {}

    with open(lc_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # 找到METRICS CLASS行和数据行
    for i, line in enumerate(lines):
        if line.startswith('## METRICS CLASS'):
            # 下一行是列名，再下一行开始是数据
            header_line = lines[i + 1].strip()
            data_line = lines[i + 2].strip()

            if data_line and not data_line.startswith('##'):
                headers = header_line.split('\t')
                values = data_line.split('\t')

                # 提取关键指标
                key_fields = [
                    'ESTIMATED_LIBRARY_SIZE',
                    'UNIQUE_READ_PAIRS',
                    'DUPLICATE_READ_PAIRS',
                    'PERCENT_DUPLICATION',
                    'MAX_LIBRARY_SIZE',
                    'MIN_LIBRARY_SIZE'
                ]

                for idx, header in enumerate(headers):
                    if header in key_fields and idx < len(values):
                        key_name = header.lower()
                        try:
                            result[key_name] = float(values[idx])
                        except ValueError:
                            result[key_name] = values[idx]
            break

    return result


def main():
    parser = argparse.ArgumentParser(description='生成全外显子质控报表')
    parser.add_argument('--sample', required=True, help='样本编号')
    parser.add_argument('--output', required=True, help='输出JSON报表路径')
    parser.add_argument('--sry-cutoff', type=int, default=10, help='SRY count cutoff for gender determination (default: 10)')
    parser.add_argument('--fastp', help='fastp_stats.json文件路径')
    parser.add_argument('--xamdst', help='xamdst.report.json文件路径')
    parser.add_argument('--mt-xamdst', help='mt.xamdst.report.json文件路径')
    parser.add_argument('--fingerprint', help='fingerprint.json文件路径')
    parser.add_argument('--metrics', help='metrics.txt文件路径')
    parser.add_argument('--hs', help='hs.txt文件路径')
    parser.add_argument('--sry', help='SRY.count.txt文件路径')
    parser.add_argument('--library-complexity', help='GATK EstimateLibraryComplexity输出文件路径')
    args = parser.parse_args()

    # 初始化报表数据
    report = {
        'sample_id': args.sample,
        'fastp': {},
        'xamdst': {},
        'mt_xamdst': {},
        'fingerprint': {},
        'metrics': {},
        'hs_metrics': {},
        'sry': {},
        'library_complexity': {}
    }

    # 解析各个输入文件
    if args.fastp:
        report['fastp'] = parse_fastp_stats(args.fastp)

    if args.xamdst:
        report['xamdst'] = parse_xamdst_report(args.xamdst)

    if args.mt_xamdst:
        report['mt_xamdst'] = parse_mt_xamdst_report(args.mt_xamdst)

    if args.fingerprint:
        report['fingerprint'] = parse_fingerprint(args.fingerprint)

    if args.metrics:
        report['metrics'] = parse_metrics(args.metrics)

    if args.hs:
        report['hs_metrics'] = parse_hs_metrics(args.hs)

    if args.sry:
        report['sry'] = parse_sry_count(args.sry, args.sry_cutoff)

    if args.library_complexity:
        report['library_complexity'] = parse_library_complexity(args.library_complexity)

    # 输出JSON报表
    with open(args.output, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)


if __name__ == '__main__':
    main()