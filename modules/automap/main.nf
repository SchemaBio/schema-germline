// AUTOMAP 模块
// 用途：ROH (Runs of Homozygosity) 分析
// 工具：AutoMap
// 输入：VCF 文件 + 参考基因组 + 捕获区域 BED
// 输出：ROH 区域检测结果 + 近亲系数估算 + 致病基因区域分析
//
// 说明：
//   - AutoMap 用于检测基因组中连续纯合区域 (ROH)
//   - ROH 分析对常染色体隐性遗传病诊断具有重要意义
//   - 可估算近亲婚配系数 (F coefficient)
//   - 支持检测单亲二倍体 (UPD) 区域
//   - 结合捕获区域 BED 文件，可分析外显子区域的 ROH

// ============================================================================
// AUTOMAP_ROH - ROH 区域检测
// ============================================================================
process AUTOMAP_ROH {
    tag "AUTOMAP on ${vcf.baseName}"
    label 'process_medium'
    label 'automap'

    input:
        path vcf                   // VCF 文件 (已注释，含样本变异)
        path vcf_tbi               // VCF 索引文件
        path fasta                 // 参考基因组 FASTA
        path fasta_fai             // 参考基因组索引 (.fai)
        path target_bed            // 捕获区域 BED 文件 (可选)
        val genome_assembly        // 基因组版本: GRCh37 或 GRCh38
        val min_roh_length         // ROH 最小长度 bp (默认 1000000 = 1Mb)
        val max_roh_length         // ROH 最大长度 bp (可选，默认不限)
        val min_snps_in_roh        // ROH 区域内最小 SNP 数 (默认 50)
        val max_gap_length         // ROH 内最大允许间隔 bp (默认 500000)
        val sample_id              // 样本 ID (可选，默认从 VCF 提取)
        val output_dir             // 输出目录 (可选)

    output:
        path "*.roh.bed", emit: roh_bed
        path "*.roh.csv", emit: roh_csv
        path "*.roh.summary.json", emit: summary_json
        path "*.roh.chromosome.txt", emit: chrom_report
        path "*.inbreeding.json", emit: inbreeding_json

    script:
    def assembly = genome_assembly ?: 'GRCh38'
    def min_length = min_roh_length ?: 1000000
    def max_length = max_roh_length ?: 100000000
    def min_snps = min_snps_in_roh ?: 50
    def max_gap = max_gap_length ?: 500000
    def sid = sample_id ?: vcf.baseName.replaceAll(/\.(vcf|vcf\.gz)$/, '')
    """
    # 运行 AutoMap ROH 分析
    # AutoMap 参数说明:
    #   --min-length: ROH 最小长度阈值，过滤短片段假阳性
    #   --max-length: ROH 最大长度阈值，用于检测不同类型的 ROH
    #   --min-snps: ROH 内最小 SNP 数，确保统计可靠性
    #   --max-gap: ROH 内允许的最大间隔，控制 ROH 连续性
    #   --target-bed: 仅分析捕获区域内的 ROH (适用于 WES)

    # 创建样本 ID 文件 (如果 VCF 包含多个样本)
    # AutoMap 需要指定分析的样本
    echo "${sid}" > samples.txt

    # 运行 AutoMap
    automap \\
        --vcf ${vcf} \\
        --fasta ${fasta} \\
        --sample ${sid} \\
        --min-length ${min_length} \\
        --max-length ${max_length} \\
        --min-snps ${min_snps} \\
        --max-gap ${max_gap} \\
        --assembly ${assembly} \\
        --output ${sid}.roh \\
        ${target_bed ? "--target-bed ${target_bed}" : ''}

    # 计算 ROH 统计信息
    # 包括：总 ROH 数量、总长度、染色体分布、近亲系数
    python3 << 'PYEOF'
    import json
    import csv
    import os

    # 读取 ROH BED 文件
    roh_regions = []
    total_roh_length = 0
    chrom_stats = {}

    roh_bed_file = "${sid}.roh.bed"
    if os.path.exists(roh_bed_file):
        with open(roh_bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    length = end - start
                    roh_regions.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'length': length,
                        'snp_count': parts[3] if len(parts) > 3 else 'NA'
                    })
                    total_roh_length += length

                    # 染色体统计
                    if chrom not in chrom_stats:
                        chrom_stats[chrom] = {'count': 0, 'total_length': 0}
                    chrom_stats[chrom]['count'] += 1
                    chrom_stats[chrom]['total_length'] += length

    # 计算近亲系数 (F coefficient)
    # F = ROH 总长度 / 基因组总长度 (近似计算)
    # 人类基因组长度约为 3,000,000,000 bp
    genome_length = 3000000000
    f_coefficient = total_roh_length / genome_length if genome_length > 0 else 0

    # 根据近亲系数判断可能的近亲关系类型
    # F < 0.01: 无近亲关系 (随机婚配)
    # F ~ 0.0625: 第一代堂表亲婚配
    # F ~ 0.125: 第二代堂表亲或叔侄婚配
    # F ~ 0.25: 兄妹或亲子婚配
    # F > 0.25: 多代近亲婚配

    def interpret_f_coefficient(f):
        if f < 0.01:
            return "No significant inbreeding detected (random mating)"
        elif f < 0.0625:
            return "Possible distant cousin mating"
        elif f < 0.09:
            return "Possible second cousin mating (F ~ 0.0156)"
        elif f < 0.125:
            return "Possible first cousin once removed or double second cousin"
        elif f < 0.1875:
            return "Possible first cousin mating (F ~ 0.0625)"
        elif f < 0.25:
            return "Possible uncle/niece or half-sibling mating (F ~ 0.125)"
        elif f < 0.375:
            return "Possible sibling or parent-offspring mating (F ~ 0.25)"
        else:
            return "Severe inbreeding, multiple generations of consanguinity"

    interpretation = interpret_f_coefficient(f_coefficient)

    # 分类 ROH 长度
    # 短 ROH (< 2 Mb): 可能是古代近亲婚配的结果
    # 中等 ROH (2-4 Mb): 可能是几代内的近亲婚配
    # 长 ROH (> 4 Mb): 可能是近几代的近亲婚配
    short_roh = [r for r in roh_regions if r['length'] < 2000000]
    medium_roh = [r for r in roh_regions if 2000000 <= r['length'] < 4000000]
    long_roh = [r for r in roh_regions if r['length'] >= 4000000]

    # 保存 ROH CSV
    with open("${sid}.roh.csv", 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['chrom', 'start', 'end', 'length', 'snp_count'])
        writer.writeheader()
        for r in sorted(roh_regions, key=lambda x: (x['chrom'], x['start'])):
            writer.writerow(r)

    # 保存染色体统计报告
    with open("${sid}.roh.chromosome.txt", 'w') as f:
        f.write("ROH Chromosome Distribution Report\\n")
        f.write("=" * 50 + "\\n\\n")
        f.write(f"Total ROH regions: {len(roh_regions)}\\n")
        f.write(f"Total ROH length: {total_roh_length} bp ({total_roh_length/1000000:.2f} Mb)\\n\\n")
        f.write("Short ROH (< 2 Mb): {} regions\\n".format(len(short_roh)))
        f.write("Medium ROH (2-4 Mb): {} regions\\n".format(len(medium_roh)))
        f.write("Long ROH (> 4 Mb): {} regions\\n\\n".format(len(long_roh)))
        f.write("Per-chromosome statistics:\\n")
        f.write("-" * 40 + "\\n")
        for chrom in sorted(chrom_stats.keys()):
            stats = chrom_stats[chrom]
            f.write(f"{chrom}: {stats['count']} ROH, {stats['total_length']/1000000:.2f} Mb\\n")

    # 保存近亲系数 JSON
    with open("${sid}.inbreeding.json", 'w') as f:
        json.dump({
            'sample_id': "${sid}",
            'f_coefficient': round(f_coefficient, 6),
            'total_roh_length': total_roh_length,
            'total_roh_length_mb': round(total_roh_length/1000000, 2),
            'num_roh_regions': len(roh_regions),
            'interpretation': interpretation,
            'roh_classification': {
                'short_roh_count': len(short_roh),
                'medium_roh_count': len(medium_roh),
                'long_roh_count': len(long_roh)
            },
            'chromosome_stats': chrom_stats
        }, f, indent=2)

    # 保存汇总 JSON
    with open("${sid}.roh.summary.json", 'w') as f:
        json.dump({
            'sample_id': "${sid}",
            'genome_assembly': "${assembly}",
            'min_roh_length': ${min_length},
            'min_snps_in_roh': ${min_snps},
            'roh_regions': roh_regions,
            'statistics': {
                'total_count': len(roh_regions),
                'total_length_bp': total_roh_length,
                'total_length_mb': round(total_roh_length/1000000, 2),
                'f_coefficient': round(f_coefficient, 6),
                'avg_roh_length': round(total_roh_length/len(roh_regions), 2) if roh_regions else 0
            }
        }, f, indent=2)
    PYEOF
    """
}


// ============================================================================
// AUTOMAP_UPD - 单亲二倍体 (UPD) 检测
// ============================================================================
process AUTOMAP_UPD {
    tag "AUTOMAP_UPD: ${sample_id}"
    label 'process_medium'
    label 'automap'

    input:
        path vcf                   // VCF 文件
        path vcf_tbi               // VCF 索引文件
        path fasta                 // 参考基因组 FASTA
        path fasta_fai             // 参考基因组索引 (.fai)
        val genome_assembly        // 基因组版本
        val sample_id              // 样本 ID
        val father_id              // 父样本 ID (可选)
        val mother_id              // 母样本 ID (可选)
        path father_vcf            // 父样本 VCF (可选)
        path father_vcf_tbi        // 父样本 VCF 索引 (可选)
        path mother_vcf            // 母样本 VCF (可选)
        path mother_vcf_tbi        // 母样本 VCF 索引 (可选)
        val output_dir             // 输出目录 (可选)

    output:
        path "*.upd.bed", emit: upd_bed
        path "*.upd.json", emit: upd_json
        path "*.upd.report.txt", emit: upd_report

    script:
    def assembly = genome_assembly ?: 'GRCh38'
    def sid = sample_id ?: vcf.baseName.replaceAll(/\.(vcf|vcf\.gz)$/, '')
    """
    # 单亲二倍体 (UPD) 检测
    # UPD 是指个体从同一亲本继承了某染色体的两个拷贝
    # 类型:
    #   - Uniparental isodisomy (UPiD): 同一染色体的两个拷贝 (导致纯合)
    #   - Uniparental heterodisomy (UPhD): 同一亲本的两个不同染色体拷贝
    # UPD 检测需要父母的 VCF 数据进行比对

    # 检查是否提供父母样本
    HAS_FATHER=false
    HAS_MOTHER=false

    if [ -n "${father_vcf}" ] && [ -f "${father_vcf}" ]; then
        HAS_FATHER=true
    fi
    if [ -n "${mother_vcf}" ] && [ -f "${mother_vcf}" ]; then
        HAS_MOTHER=true
    fi

    # 运行 AutoMap UPD 检测
    # 如果有父母数据，可以进行更精确的 UPD 检测
    # 如果没有父母数据，仅基于 ROH 分析推断可能的 UPD 区域

    if [ "\$HAS_FATHER" = true ] && [ "\$HAS_MOTHER" = true ]; then
        # Trio 模式：有父母数据，精确 UPD 检测
        automap-upd \\
            --vcf ${vcf} \\
            --father-vcf ${father_vcf} \\
            --mother-vcf ${mother_vcf} \\
            --sample ${sid} \\
            --father-id ${father_id} \\
            --mother-id ${mother_id} \\
            --fasta ${fasta} \\
            --assembly ${assembly} \\
            --output ${sid}.upd
    else
        # 单样本模式：仅基于 ROH 推断可能的 UPD 区域
        # ROH 长染色体 (> 10 Mb) 可能指示 UPD
        automap-upd \\
            --vcf ${vcf} \\
            --sample ${sid} \\
            --fasta ${fasta} \\
            --assembly ${assembly} \\
            --infer-mode \\
            --output ${sid}.upd
    fi

    # 生成 UPD 报告
    python3 << 'PYEOF'
    import json
    import os

    # 读取 UPD 结果
    upd_regions = []
    upd_bed_file = "${sid}.upd.bed"

    if os.path.exists(upd_bed_file):
        with open(upd_bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    upd_regions.append({
                        'chrom': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'upd_type': parts[3] if len(parts) > 3 else 'unknown'
                    })

    # 生成报告
    with open("${sid}.upd.report.txt", 'w') as f:
        f.write("UPD Detection Report for ${sid}\\n")
        f.write("=" * 50 + "\\n\\n")

        if not upd_regions:
            f.write("No UPD regions detected.\\n")
            f.write("All chromosomes show biparental inheritance pattern.\\n")
        else:
            f.write(f"Detected {len(upd_regions)} potential UPD regions:\\n\\n")
            for region in upd_regions:
                f.write(f"Chromosome {region['chrom']}:\\n")
                f.write(f"  Region: {region['start']}-{region['end']} ({(region['end']-region['start'])/1000000:.2f} Mb)\\n")
                f.write(f"  Type: {region['upd_type']}\\n\\n")

            f.write("\\nClinical significance:\\n")
            f.write("-" * 30 + "\\n")
            f.write("UPD regions may lead to:\\n")
            f.write("  1. Homozygosity for recessive disease variants\\n")
            f.write("  2. Imprinting disorders (e.g., Prader-Willi, Angelman)\\n")
            f.write("  3. Growth and developmental abnormalities\\n")
            f.write("\\nRecommended follow-up:\\n")
            f.write("  - Verify with parental samples if available\\n")
            f.write("  - Check for known imprinting disorders on affected chromosomes\\n")
            f.write("  - Analyze genes within UPD regions for recessive variants\\n")

    # 保存 JSON
    with open("${sid}.upd.json", 'w') as f:
        json.dump({
            'sample_id': "${sid}",
            'genome_assembly': "${assembly}",
            'has_parental_data': ${father_vcf ? 'true' : 'false'} and ${mother_vcf ? 'true' : 'false'},
            'upd_regions': upd_regions,
            'total_upd_regions': len(upd_regions)
        }, f, indent=2)
    PYEOF
    """
}


// ============================================================================
// AUTOMAP_GENE_OVERLAY - ROH 区域致病基因分析
// ============================================================================
process AUTOMAP_GENE_OVERLAY {
    tag "AUTOMAP_GENE: ${sample_id}"
    label 'process_low'
    label 'automap'

    input:
        path roh_bed               // ROH BED 文件
        path gene_bed              // 基因 BED 文件 (基因区域定义)
        path disease_genes         // 致病基因列表文件 (可选，JSON/CSV 格式)
        val sample_id              // 样本 ID
        val min_gene_overlap       // 最小基因重叠比例 (默认 0.5)
        val output_dir             // 输出目录 (可选)

    output:
        path "*.roh.genes.csv", emit: genes_csv
        path "*.roh.disease_genes.csv", emit: disease_genes_csv
        path "*.roh.genes.json", emit: genes_json
        path "*.roh.genes.report.txt", emit: genes_report

    script:
    def sid = sample_id ?: "sample"
    def min_overlap = min_gene_overlap ?: 0.5
    """
    # ROH 区域致病基因分析
    # 将 ROH 区域与基因区域进行叠加分析
    # 找出完全或部分位于 ROH 区域内的基因
    # 特别关注常染色体隐性遗传病的致病基因

    # 使用 bedtools 进行区域叠加分析
    # bedtools intersect 计算 ROH 与基因的重叠

    # 找出位于 ROH 区域内的所有基因
    bedtools intersect \\
        -a ${gene_bed} \\
        -b ${roh_bed} \\
        -wa -wb \\
        -F ${min_overlap} \\
        > ${sid}.roh.genes.intersect.bed

    # 分析重叠基因
    python3 << 'PYEOF'
    import json
    import csv
    import os

    # 读取基因重叠结果
    genes_in_roh = []

    intersect_file = "${sid}.roh.genes.intersect.bed"
    if os.path.exists(intersect_file):
        with open(intersect_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    # 基因信息 (from gene_bed)
                    gene_chrom = parts[0]
                    gene_start = int(parts[1])
                    gene_end = int(parts[2])
                    gene_name = parts[3] if len(parts) > 3 else 'unknown'
                    gene_length = gene_end - gene_start

                    # ROH 信息 (from roh_bed)
                    roh_chrom = parts[4]
                    roh_start = int(parts[5])
                    roh_end = int(parts[6])

                    # 计算重叠
                    overlap_start = max(gene_start, roh_start)
                    overlap_end = min(gene_end, roh_end)
                    overlap_length = overlap_end - overlap_start
                    overlap_fraction = overlap_length / gene_length if gene_length > 0 else 0

                    genes_in_roh.append({
                        'gene_name': gene_name,
                        'chrom': gene_chrom,
                        'gene_start': gene_start,
                        'gene_end': gene_end,
                        'roh_start': roh_start,
                        'roh_end': roh_end,
                        'overlap_length': overlap_length,
                        'overlap_fraction': round(overlap_fraction, 3),
                        'gene_length': gene_length
                    })

    # 读取致病基因列表 (如果提供)
    disease_genes_list = []
    disease_genes_in_roh = []

    disease_file = "${disease_genes}"
    if disease_file and os.path.exists(disease_file):
        try:
            with open(disease_file, 'r') as f:
                if disease_file.endswith('.json'):
                    data = json.load(f)
                    disease_genes_list = data.get('genes', [])
                else:
                    # CSV 格式
                    reader = csv.DictReader(f)
                    for row in reader:
                        disease_genes_list.append(row.get('gene', row.get('gene_name', '')))
        except Exception as e:
            print(f"Warning: Could not read disease genes file: {e}")

    # 找出位于 ROH 的致病基因
    for gene in genes_in_roh:
        if gene['gene_name'] in disease_genes_list:
            disease_genes_in_roh.append(gene)

    # 保存结果
    # 所有基因 CSV
    with open("${sid}.roh.genes.csv", 'w', newline='') as f:
        if genes_in_roh:
            writer = csv.DictWriter(f, fieldnames=genes_in_roh[0].keys())
            writer.writeheader()
            writer.writerows(genes_in_roh)

    # 致病基因 CSV
    with open("${sid}.roh.disease_genes.csv", 'w', newline='') as f:
        if disease_genes_in_roh:
            writer = csv.DictWriter(f, fieldnames=disease_genes_in_roh[0].keys())
            writer.writeheader()
            writer.writerows(disease_genes_in_roh)

    # 保存 JSON
    with open("${sid}.roh.genes.json", 'w') as f:
        json.dump({
            'sample_id': "${sid}",
            'min_overlap_fraction': ${min_overlap},
            'total_genes_in_roh': len(genes_in_roh),
            'genes_in_roh': genes_in_roh,
            'disease_genes_in_roh': disease_genes_in_roh,
            'num_disease_genes_in_roh': len(disease_genes_in_roh)
        }, f, indent=2)

    # 生成报告
    with open("${sid}.roh.genes.report.txt", 'w') as f:
        f.write("ROH Gene Analysis Report for ${sid}\\n")
        f.write("=" * 50 + "\\n\\n")
        f.write(f"Total genes overlapping ROH regions: {len(genes_in_roh)}\\n")
        f.write(f"Disease genes in ROH regions: {len(disease_genes_in_roh)}\\n\\n")

        if disease_genes_in_roh:
            f.write("Disease genes located in ROH regions:\\n")
            f.write("-" * 40 + "\\n")
            for gene in disease_genes_in_roh:
                f.write(f"Gene: {gene['gene_name']}\\n")
                f.write(f"  Chromosome: {gene['chrom']}\\n")
                f.write(f"  Gene position: {gene['gene_start']}-{gene['gene_end']}\\n")
                f.write(f"  ROH region: {gene['roh_start']}-{gene['roh_end']}\\n")
                f.write(f"  Overlap: {gene['overlap_length']} bp ({gene['overlap_fraction']*100:.1f}%)\\n\\n")

            f.write("\\nClinical implications:\\n")
            f.write("-" * 30 + "\\n")
            f.write("Genes in ROH regions are homozygous.\\n")
            f.write("Recessive disease variants in these genes may be expressed.\\n")
            f.write("Consider analyzing VCF for homozygous variants in these genes.\\n")
        else:
            f.write("No known disease genes detected in ROH regions.\\n")
    PYEOF

    # 清理临时文件
    rm -f ${sid}.roh.genes.intersect.bed
    """
}


// ============================================================================
// AUTOMAP_VISUALIZE - ROH 可视化
// ============================================================================
process AUTOMAP_VISUALIZE {
    tag "AUTOMAP_VIS: ${sample_id}"
    label 'process_low'
    label 'automap'

    input:
        path roh_bed               // ROH BED 文件
        path roh_summary           // ROH 摘要 JSON 文件
        val sample_id              // 样本 ID
        val genome_assembly        // 基因组版本
        val output_dir             // 输出目录 (可选)

    output:
        path "*.roh.plot.png", emit: roh_plot
        path "*.roh.plot.html", emit: roh_html

    script:
    def sid = sample_id ?: "sample"
    def assembly = genome_assembly ?: 'GRCh38'
    """
    # ROH 可视化
    # 生成染色体级别的 ROH 分布图

    # 使用 Python matplotlib 绘制 ROH 分布图
    python3 << 'PYEOF'
    import json
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import numpy as np
    import os

    # 染色体长度参考 (GRCh38)
    chrom_lengths = {
        '1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555,
        '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
        '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
        '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
        '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
        '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415
    }

    # GRCh37 染色体长度
    chrom_lengths_hg19 = {
        '1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276,
        '5': 180915260, '6': 171115067, '7': 159138663, '8': 146364022,
        '9': 141213431, '10': 135534747, '11': 135006516, '12': 133851895,
        '13': 115169878, '14': 107349540, '15': 102531392, '16': 90354753,
        '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
        '21': 48129895, '22': 51304566, 'X': 155270560, 'Y': 59373566
    }

    # 选择参考
    if "${assembly}" == "GRCh37":
        ref_lengths = chrom_lengths_hg19
    else:
        ref_lengths = chrom_lengths

    # 读取 ROH 数据
    roh_regions = []
    with open("${roh_bed}", 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0].replace('chr', '')
                start = int(parts[1])
                end = int(parts[2])
                roh_regions.append({'chrom': chrom, 'start': start, 'end': end})

    # 读取摘要信息
    summary = {}
    if os.path.exists("${roh_summary}"):
        with open("${roh_summary}", 'r') as f:
            summary = json.load(f)

    # 创建染色体分布图
    fig, ax = plt.subplots(figsize=(16, 8))

    # 常染色体列表
    autosomes = [str(i) for i in range(1, 23)]

    # 绘制染色体背景
    for i, chrom in enumerate(autosomes):
        y = i
        chrom_len = ref_lengths.get(chrom, 100000000)

        # 染色体背景
        ax.add_patch(patches.Rectangle((0, y-0.4), chrom_len/1e6, 0.8,
                                        facecolor='lightgray', edgecolor='black'))

        # 标注染色体名称
        ax.text(-5, y, chrom, ha='right', va='center', fontsize=10)

    # 绘制 ROH 区域
    for roh in roh_regions:
        chrom = roh['chrom'].replace('chr', '')
        if chrom in autosomes:
            y = int(chrom) - 1
            start_mb = roh['start'] / 1e6
            end_mb = roh['end'] / 1e6
            length_mb = end_mb - start_mb

            # ROH 区域用红色标记
            ax.add_patch(patches.Rectangle((start_mb, y-0.35), length_mb, 0.7,
                                            facecolor='red', edgecolor='darkred',
                                            alpha=0.7))

    # 设置坐标轴
    ax.set_xlim(-10, max([ref_lengths.get(c, 100000000)/1e6 for c in autosomes]) + 10)
    ax.set_ylim(-1, 23)
    ax.set_xlabel('Position (Mb)', fontsize=12)
    ax.set_ylabel('Chromosome', fontsize=12)
    ax.set_title('ROH Distribution for ${sid}\\n'
                 f'Total ROH: {summary.get("statistics", {}).get("total_length_mb", 0):.2f} Mb, '
                 f'F coefficient: {summary.get("statistics", {}).get("f_coefficient", 0):.4f}',
                 fontsize=14)

    # 添加网格
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig("${sid}.roh.plot.png", dpi=150)
    plt.close()

    # 创建 HTML 可视化报告
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>ROH Analysis Report - ${sid}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1 {{ color: #2c3e50; }}
            h2 {{ color: #34495e; border-bottom: 1px solid #bdc3c7; }}
            .summary {{ background: #ecf0f1; padding: 15px; border-radius: 5px; }}
            img {{ max-width: 100%; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #bdc3c7; padding: 8px; text-align: left; }}
            th {{ background: #34495e; color: white; }}
        </style>
    </head>
    <body>
        <h1>ROH Analysis Report</h1>
        <p>Sample: ${sid}</p>
        <p>Assembly: ${assembly}</p>

        <h2>Summary</h2>
        <div class="summary">
            <p>Total ROH regions: {len(roh_regions)}</p>
            <p>Total ROH length: {summary.get('statistics', {}).get('total_length_mb', 0):.2f} Mb</p>
            <p>Inbreeding coefficient (F): {summary.get('statistics', {}).get('f_coefficient', 0):.4f}</p>
        </div>

        <h2>ROH Distribution Plot</h2>
        <img src="${sid}.roh.plot.png" alt="ROH Distribution">

        <h2>ROH Regions</h2>
        <table>
            <tr><th>Chromosome</th><th>Start (Mb)</th><th>End (Mb)</th><th>Length (Mb)</th></tr>
    """

    for roh in sorted(roh_regions, key=lambda x: (x['chrom'], x['start'])):
        html_content += f"""
            <tr>
                <td>{roh['chrom']}</td>
                <td>{roh['start']/1e6:.2f}</td>
                <td>{roh['end']/1e6:.2f}</td>
                <td>{(roh['end']-roh['start'])/1e6:.2f}</td>
            </tr>
        """

    html_content += """
        </table>
    </body>
    </html>
    """

    with open("${sid}.roh.plot.html", 'w') as f:
        f.write(html_content)
    PYEOF
    """
}