#!/usr/bin/env python3
"""
CNV Segmentation Script

从外显子水平的覆盖度数据计算染色体片段水平的 CNV，
并使用 cytoband 注释生成 HGVS 格式描述。

输入格式 (支持 CNVkit .cnr 或类似格式):
    chromosome  start   end     gene    depth   weight  log2
    chr1        10000   10500   BRCA1   50      1       0.0
    chr1        11000   11500   TP53    48      1       -0.05

输出格式:
    chromosome  start   end     log2    copies  cnv_type     cytoband            hgvs

算法: 基于统计检验的 Circular Binary Segmentation (CBS) 简化实现

Usage:
    python cnv_segment.py -i sample.cnr -o sample.seg.tsv
    python cnv_segment.py -i sample.cnr -o sample.seg.tsv --cytoband cytoband.txt
    python cnv_segment.py -i sample.cnr -o sample.seg.tsv --cytoband cytoband.txt --assembly GRCh38
"""

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import numpy as np
from scipy import stats


# 参考基因组 ID 映射
REFSEQ_TO_CHROM = {
    'GRCh37': {
        'chr1': 'NC_000001.10', 'chr2': 'NC_000002.11', 'chr3': 'NC_000003.12',
        'chr4': 'NC_000004.12', 'chr5': 'NC_000005.10', 'chr6': 'NC_000006.12',
        'chr7': 'NC_000007.14', 'chr8': 'NC_000008.11', 'chr9': 'NC_000009.12',
        'chr10': 'NC_000010.11', 'chr11': 'NC_000011.10', 'chr12': 'NC_000012.12',
        'chr13': 'NC_000013.11', 'chr14': 'NC_000014.9', 'chr15': 'NC_000015.10',
        'chr16': 'NC_000016.11', 'chr17': 'NC_000017.11', 'chr18': 'NC_000018.10',
        'chr19': 'NC_000019.10', 'chr20': 'NC_000020.11', 'chr21': 'NC_000021.9',
        'chr22': 'NC_000022.10', 'chrX': 'NC_000023.11', 'chrY': 'NC_000024.10'
    },
    'GRCh38': {
        'chr1': 'NC_000001.11', 'chr2': 'NC_000002.12', 'chr3': 'NC_000003.12',
        'chr4': 'NC_000004.12', 'chr5': 'NC_000005.10', 'chr6': 'NC_000006.12',
        'chr7': 'NC_000007.14', 'chr8': 'NC_000008.11', 'chr9': 'NC_000009.12',
        'chr10': 'NC_000010.11', 'chr11': 'NC_000011.10', 'chr12': 'NC_000012.12',
        'chr13': 'NC_000013.11', 'chr14': 'NC_000014.9', 'chr15': 'NC_000015.10',
        'chr16': 'NC_000016.11', 'chr17': 'NC_000017.11', 'chr18': 'NC_000018.10',
        'chr19': 'NC_000019.10', 'chr20': 'NC_000020.11', 'chr21': 'NC_000021.9',
        'chr22': 'NC_000022.11', 'chrX': 'NC_000023.11', 'chrY': 'NC_000024.10'
    }
}


@dataclass
class ExonRecord:
    """外显子记录"""
    chrom: str
    start: int
    end: int
    gene: Optional[str]
    depth: float
    weight: float
    log2: float


@dataclass
class Cytoband:
    """细胞色素带"""
    chrom: str
    start: int
    end: int
    name: str
    stain: str


@dataclass
class Segment:
    """分段结果"""
    chrom: str
    start: int
    end: int
    log2: float
    copies: int
    cnv_type: str
    probes: int
    genes: list[str]
    cytobands: list[str] = None
    hgvs: str = None


def parse_cnr(filepath: Path) -> list[ExonRecord]:
    """解析 CNVkit .cnr 或类似格式文件"""
    records = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 7:
                continue

            try:
                record = ExonRecord(
                    chrom=parts[0],
                    start=int(parts[1]),
                    end=int(parts[2]),
                    gene=parts[3] if parts[3] else None,
                    depth=float(parts[4]),
                    weight=float(parts[5]) if len(parts) > 5 else 1.0,
                    log2=float(parts[6]) if len(parts) > 6 else 0.0
                )
                records.append(record)
            except (ValueError, IndexError):
                continue

    return records


def parse_cytoband(filepath: Path) -> list[Cytoband]:
    """解析 cytoband 文件 (UCSC 格式)"""
    cytobands = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 5:
                continue

            try:
                band = Cytoband(
                    chrom=parts[0],
                    start=int(parts[1]),
                    end=int(parts[2]),
                    name=parts[3],
                    stain=parts[4]
                )
                cytobands.append(band)
            except (ValueError, IndexError):
                continue

    return cytobands


def build_cytoband_index(cytobands: list[Cytoband]) -> dict[str, list[Cytoband]]:
    """按染色体建立 cytoband 索引"""
    index = {}
    for band in cytobands:
        if band.chrom not in index:
            index[band.chrom] = []
        index[band.chrom].append(band)
    return index


def get_cytobands_for_region(chrom: str, start: int, end: int,
                               cytoband_index: dict[str, list[Cytoband]]) -> list[str]:
    """获取一个区域跨越的所有 cytoband 名称"""
    if chrom not in cytoband_index:
        return []

    bands = cytoband_index[chrom]
    result = []

    for band in bands:
        # 检查重叠
        if band.end <= start or band.start >= end:
            continue
        result.append(band.name)

    return result


def generate_hgvs(chrom: str, start: int, end: int, cnv_type: str,
                   assembly: str = 'GRCh38') -> str:
    """
    生成 HGVS 格式描述

    HGVS 格式规则:
    - 缺失: chr1:g.100000_200000del
    - 重复: chr1:g.100000_200000dup

    Args:
        chrom: 染色体名称 (如 chr1)
        start: 起始位置 (1-based)
        end: 结束位置 (1-based, inclusive)
        cnv_type: CNV 类型 (del, dup, gain, loss, amplification)
        assembly: 基因组版本 (当前未使用，保留兼容性)

    Returns:
        HGVS 格式描述字符串
    """
    # 确定变异类型
    if cnv_type in ('loss', 'heterozygous_loss', 'neutral'):
        if cnv_type == 'neutral':
            return '.'
        variant_type = 'del'
    elif cnv_type in ('gain', 'amplification'):
        variant_type = 'dup'
    else:
        variant_type = cnv_type

    # 使用 chr 前缀格式
    return f"{chrom}:g.{start}_{end}{variant_type}"


def generate_cytoband_description(cytoband_names: list[str], cnv_type: str) -> str:
    """
    生成基于 cytoband 的临床描述

    例如: 1q21.1q21.2 duplication
          22q11.2 deletion

    Args:
        cytoband_names: 跨越的 cytoband 列表
        cnv_type: CNV 类型

    Returns:
        临床描述字符串
    """
    if not cytoband_names:
        return '.'

    # 解析 cytoband 名称，提取臂和区
    # cytoband 格式: p36.33, q21.1, etc.
    arm_bands = {}
    for name in cytoband_names:
        if name.startswith('p'):
            arm = 'p'
        elif name.startswith('q'):
            arm = 'q'
        else:
            arm = ''

        # 提取区号 (如 21.1 -> 21.1)
        region = name[1:] if arm else name
        if arm not in arm_bands:
            arm_bands[arm] = []
        arm_bands[arm].append(region)

    # 构建描述
    descriptions = []
    for arm in ['p', 'q']:
        if arm in arm_bands and arm_bands[arm]:
            bands = sorted(arm_bands[arm], key=lambda x: (
                float(x.split('.')[0]) if '.' in x and x.split('.')[0].isdigit() else 0,
                float(x.split('.')[1]) if '.' in x and len(x.split('.')) > 1 else 0
            ))

            if len(bands) == 1:
                desc = f"{arm}{bands[0]}"
            else:
                # 找到起始和结束的区
                start_band = bands[0]
                end_band = bands[-1]
                # 提取共同前缀
                common_prefix = ''
                for i, c in enumerate(start_band):
                    if i < len(end_band) and start_band[i] == end_band[i]:
                        common_prefix += c
                    else:
                        break
                # 移除末尾的不完整数字
                while common_prefix and (common_prefix[-1].isdigit() and
                       (len(common_prefix) == len(start_band) or
                        common_prefix[-1] != start_band[len(common_prefix)])):
                    common_prefix = common_prefix[:-1]

                desc = f"{arm}{common_prefix}"

                # 如果结束区有额外的后缀部分，添加
                end_suffix = end_band[len(common_prefix):] if len(end_band) > len(common_prefix) else ''
                if end_suffix:
                    desc += end_suffix

            descriptions.append(desc)

    if not descriptions:
        return cytoband_names[0] if cytoband_names else '.'

    # 添加变异类型
    if cnv_type in ('loss', 'heterozygous_loss'):
        cnv_word = 'deletion'
    elif cnv_type in ('gain', 'amplification'):
        cnv_word = 'duplication'
    elif cnv_type == 'neutral':
        return '.'
    else:
        cnv_word = cnv_type

    if len(descriptions) == 1:
        return f"{descriptions[0]} {cnv_word}"
    else:
        return f"{'.'.join(descriptions)} {cnv_word}"


def group_by_chromosome(records: list[ExonRecord]) -> dict[str, list[ExonRecord]]:
    """按染色体分组"""
    groups = {}
    for record in records:
        if record.chrom not in groups:
            groups[record.chrom] = []
        groups[record.chrom].append(record)
    return groups


def cbs_segmentation(log2_values: np.ndarray,
                      weights: Optional[np.ndarray] = None,
                      alpha: float = 0.01,
                      min_probes: int = 3) -> list[tuple[int, int, float]]:
    """
    Circular Binary Segmentation (CBS) 简化实现

    Args:
        log2_values: log2 覆盖度值数组
        weights: 权重数组
        alpha: 显著性水平
        min_probes: 最小探针数

    Returns:
        分段列表 [(start_idx, end_idx, mean_log2), ...]
    """
    n = len(log2_values)
    if n < min_probes * 2:
        return [(0, n - 1, np.mean(log2_values))]

    if weights is None:
        weights = np.ones(n)

    segments = [(0, n - 1)]

    while True:
        best_gain = 0
        best_breakpoint = -1
        best_segment_idx = -1

        for seg_idx, (start, end) in enumerate(segments):
            seg_length = end - start + 1
            if seg_length < min_probes * 2:
                continue

            seg_values = log2_values[start:end + 1]
            seg_weights = weights[start:end + 1]

            for bp in range(start + min_probes, end - min_probes + 1):
                left_values = log2_values[start:bp]
                left_weights = weights[start:bp]
                right_values = log2_values[bp:end + 1]
                right_weights = weights[bp:end + 1]

                if len(left_values) < 2 or len(right_values) < 2:
                    continue

                left_mean = np.average(left_values, weights=left_weights)
                right_mean = np.average(right_values, weights=right_weights)

                try:
                    _, p_value = stats.ttest_ind(left_values, right_values, equal_var=False)
                except Exception:
                    continue

                gain = abs(left_mean - right_mean)

                if p_value < alpha and gain > best_gain:
                    best_gain = gain
                    best_breakpoint = bp
                    best_segment_idx = seg_idx

        if best_breakpoint == -1:
            break

        start, end = segments[best_segment_idx]
        segments.pop(best_segment_idx)
        segments.append((start, best_breakpoint - 1))
        segments.append((best_breakpoint, end))

    final_segments = []
    for start, end in segments:
        seg_values = log2_values[start:end + 1]
        seg_weights = weights[start:end + 1]
        mean_log2 = np.average(seg_values, weights=seg_weights)
        final_segments.append((start, end, mean_log2))

    return final_segments


def segment_chromosome(records: list[ExonRecord],
                       alpha: float = 0.01,
                       min_probes: int = 3) -> list[Segment]:
    """对单个染色体进行分段"""
    if len(records) < min_probes:
        log2_values = np.array([r.log2 for r in records])
        return [Segment(
            chrom=records[0].chrom,
            start=records[0].start,
            end=records[-1].end,
            log2=float(np.mean(log2_values)),
            copies=2,
            cnv_type='neutral',
            probes=len(records),
            genes=list(set(r.gene for r in records if r.gene))
        )]

    log2_values = np.array([r.log2 for r in records])
    weights = np.array([r.weight for r in records])

    valid_mask = np.isfinite(log2_values)
    if not np.all(valid_mask):
        log2_values = log2_values[valid_mask]
        weights = weights[valid_mask]
        records = [records[i] for i in np.where(valid_mask)[0]]

    if len(log2_values) < min_probes * 2:
        mean_log2 = float(np.mean(log2_values))
        return [Segment(
            chrom=records[0].chrom,
            start=records[0].start,
            end=records[-1].end,
            log2=mean_log2,
            copies=2,
            cnv_type='neutral',
            probes=len(records),
            genes=list(set(r.gene for r in records if r.gene))
        )]

    segments_idx = cbs_segmentation(log2_values, weights=weights, alpha=alpha, min_probes=min_probes)

    results = []
    for start_idx, end_idx, mean_log2 in segments_idx:
        seg_records = records[start_idx:end_idx + 1]

        copies = round(2 * (2 ** mean_log2))
        copies = max(1, min(10, copies))

        if copies == 2:
            if abs(mean_log2) < 0.2:
                cnv_type = 'neutral'
            elif mean_log2 > 0.2:
                cnv_type = 'gain'
            else:
                cnv_type = 'loss'
        elif copies > 2:
            cnv_type = 'amplification' if copies >= 5 else 'gain'
        else:
            cnv_type = 'heterozygous_loss'

        genes = list(set(r.gene for r in seg_records if r.gene))

        results.append(Segment(
            chrom=seg_records[0].chrom,
            start=seg_records[0].start,
            end=seg_records[-1].end,
            log2=mean_log2,
            copies=copies,
            cnv_type=cnv_type,
            probes=len(seg_records),
            genes=genes
        ))

    return results


def annotate_segments(segments: list[Segment],
                      cytoband_index: Optional[dict[str, list[Cytoband]]],
                      assembly: str = 'GRCh38'):
    """为分段添加 cytoband 注释和 HGVS"""
    for seg in segments:
        if cytoband_index:
            seg.cytobands = get_cytobands_for_region(seg.chrom, seg.start, seg.end, cytoband_index)
        else:
            seg.cytobands = []

        seg.cytoband_desc = generate_cytoband_description(seg.cytobands, seg.cnv_type)
        seg.hgvs = generate_hgvs(seg.chrom, seg.start, seg.end, seg.cnv_type, assembly)


def write_segments(segments: list[Segment], output_file: Path):
    """写入分段结果"""
    with open(output_file, 'w') as f:
        f.write('# CNV Segmentation Results\n')
        f.write('# Generated by cnv_segment.py\n')
        f.write('# \n')
        f.write('chromosome\tstart\tend\tlog2\tcopies\tcnv_type\tprobes\tgenes\tcytoband\thgvs\n')

        for seg in segments:
            genes_str = ','.join(seg.genes) if seg.genes else '.'
            cytoband_str = ','.join(seg.cytobands) if seg.cytobands else '.'
            hgvs_str = seg.hgvs if seg.hgvs else '.'

            f.write(f'{seg.chrom}\t{seg.start}\t{seg.end}\t'
                    f'{seg.log2:.4f}\t{seg.copies}\t{seg.cnv_type}\t'
                    f'{seg.probes}\t{genes_str}\t{cytoband_str}\t{hgvs_str}\n')


def write_bed(segments: list[Segment], output_file: Path):
    """写入 BED 格式 (用于 IGV 等可视化)"""
    with open(output_file, 'w') as f:
        f.write('#track name="CNV" description="CNV segments" visibility=2 itemRgb=On\n')

        for seg in segments:
            cnv_color = {
                'neutral': '0,128,0',
                'gain': '255,128,0',
                'amplification': '255,0,0',
                'loss': '0,0,255',
                'heterozygous_loss': '0,0,200'
            }.get(seg.cnv_type, '128,128,128')

            name = f"{seg.cnv_type}_{seg.copies}"
            score = min(1000, int(abs(seg.log2) * 100))

            f.write(f'{seg.chrom}\t{seg.start}\t{seg.end}\t'
                    f'{name}\t{score}\t.\t{cnv_color}\n')


def main():
    parser = argparse.ArgumentParser(
        description='CNV Segmentation: 从外显子水平覆盖度数据计算染色体片段水平 CNV'
    )
    parser.add_argument('-i', '--input', required=True, type=Path,
                        help='输入文件 (CNVkit .cnr 或类似格式)')
    parser.add_argument('-o', '--output', required=True, type=Path,
                        help='输出分段结果 (TSV 格式)')
    parser.add_argument('--cytoband', type=Path,
                        help='cytoband 文件 (UCSC 格式)')
    parser.add_argument('--assembly', type=str, default='GRCh38',
                        choices=['GRCh37', 'GRCh38'],
                        help='参考基因组版本 (default: GRCh38)')
    parser.add_argument('--alpha', type=float, default=0.01,
                        help='统计显著性水平 (default: 0.01)')
    parser.add_argument('--min-probes', type=int, default=3,
                        help='最小探针数 (default: 3)')
    parser.add_argument('--bed', action='store_true',
                        help='同时输出 BED 格式文件')
    parser.add_argument('--json', action='store_true',
                        help='同时输出 JSON 格式文件')

    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # 解析 cytoband 文件
    cytoband_index = None
    if args.cytoband and args.cytoband.exists():
        print(f"Parsing cytoband file: {args.cytoband}")
        cytobands = parse_cytoband(args.cytoband)
        cytoband_index = build_cytoband_index(cytobands)
        print(f"Loaded {len(cytobands)} cytobands")
    elif args.cytoband:
        print(f"Warning: Cytoband file not found: {args.cytoband}", file=sys.stderr)

    print(f"Parsing input file: {args.input}")
    records = parse_cnr(args.input)
    print(f"Loaded {len(records)} exon records")

    if not records:
        print("Error: No valid records found in input file", file=sys.stderr)
        sys.exit(1)

    print(f"Segmenting chromosomes with alpha={args.alpha}, min_probes={args.min_probes}")
    chrom_groups = group_by_chromosome(records)

    all_segments = []
    for chrom in sorted(chrom_groups.keys(), key=lambda x: (
        x.startswith('chr') and x[3:].isdigit() or not x[3:].isdigit(),
        x[3:] if x.startswith('chr') and x[3:].isdigit() else x
    )):
        chrom_records = chrom_groups[chrom]
        print(f"  Processing {chrom}: {len(chrom_records)} probes")
        chrom_segments = segment_chromosome(
            chrom_records,
            alpha=args.alpha,
            min_probes=args.min_probes
        )
        all_segments.extend(chrom_segments)

    # 添加 cytoband 注释和 HGVS
    print(f"Annotating segments with cytoband and HGVS (assembly: {args.assembly})")
    annotate_segments(all_segments, cytoband_index, args.assembly)

    # 过滤 neutral 片段（可选）
    non_neutral = [s for s in all_segments if s.cnv_type != 'neutral']
    print(f"Found {len(all_segments)} segments ({len(non_neutral)} non-neutral)")

    # 输出结果
    write_segments(all_segments, args.output)
    print(f"Written segmentation results to {args.output}")

    if args.bed:
        bed_file = args.output.with_suffix('.bed')
        write_bed(all_segments, bed_file)
        print(f"Written BED file to {bed_file}")

    if args.json:
        import json
        json_file = args.output.with_suffix('.json')
        json_data = [
            {
                'chromosome': seg.chrom,
                'start': seg.start,
                'end': seg.end,
                'log2': round(seg.log2, 4),
                'copies': seg.copies,
                'cnv_type': seg.cnv_type,
                'probes': seg.probes,
                'genes': seg.genes,
                'cytobands': seg.cytobands,
                'cytoband_description': seg.cytoband_desc,
                'hgvs': seg.hgvs
            }
            for seg in all_segments
        ]
        with open(json_file, 'w') as f:
            json.dump(json_data, f, indent=2)
        print(f"Written JSON file to {json_file}")


if __name__ == '__main__':
    main()
