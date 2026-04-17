// WHATSHAP 模块
// 用途：单倍型定相 (Haplotype Phasing)
// 工具：WhatsHap
// 输出：定相后的 VCF 文件
//
// 说明：
//   - WhatsHap 是基于 reads 的单倍型定相工具
//   - 利用测序 reads 中的变异信息进行定相
//   - 支持 BAM 输入格式
//   - 定相结果可用于确定变异的单倍型来源
//   - 常用于家系分析、亲缘关系鉴定等场景

// ============================================================================
// WHATSHAP_PHASE - 单倍型定相
// ============================================================================
process WHATSHAP_PHASE {
    tag "WHATSHAP on ${alignment.baseName}"
    label 'process_medium'
    label 'whatshap'

    input:
        path vcf                   // 输入 VCF 文件 (未定相的变异)
        path vcf_tbi               // VCF 索引
        path alignment             // BAM 比对文件
        path alignment_index       // 比对文件索引
        path fasta                 // 参考基因组 FASTA
        path fasta_fai             // 参考基因组索引 (.fai)
        val sample_id              // 样本标识符
        val chromosomes            // 要定相的染色体 (逗号分隔, 默认全部)
        val ignore_read_groups     // 是否忽略 read groups (boolean, 默认 true)
        val reference_confidence   // 参考置信度阈值 (默认 20)
        val output_dir             // 输出目录

    output:
        path "${sample_id}.phased.vcf.gz", emit: vcf
        path "${sample_id}.phased.vcf.gz.tbi", emit: vcf_tbi
        path "${sample_id}.phasing.json", emit: json, optional: true

    when:
    output_dir != 'NO_OUTPUT'

    script:
    def chrom_param = chromosomes ? "--chromosome ${chromosomes.replaceAll(/,/, ' --chromosome ')}" : ''
    def ignore_rg = ignore_read_groups ? '--ignore-read-groups' : ''
    def ref_conf = reference_confidence ?: 20
    """
    # 运行 WhatsHap 定相
    whatshap phase \\
        --reference ${fasta} \\
        --output "${sample_id}.phased.vcf" \\
        ${ignore_rg} \\
        --reference-confidence ${ref_conf} \\
        ${chrom_param} \\
        ${vcf} \\
        ${alignment}

    # 压缩 VCF 并创建索引
    bgzip -c "${sample_id}.phased.vcf" > "${sample_id}.phased.vcf.gz"
    tabix -p vcf "${sample_id}.phased.vcf.gz"

    # 清理原始 VCF
    rm -f "${sample_id}.phased.vcf"

    # 生成定相统计信息
    whatshap stats --json "${sample_id}.phasing.json" "${sample_id}.phased.vcf.gz" || true
    """
}

// ============================================================================
// WHATSHAP_HAPLOTAG - 单倍型标记
// ============================================================================
process WHATSHAP_HAPLOTAG {
    tag "HAPLOTAG on ${alignment.baseName}"
    label 'process_medium'
    label 'whatshap'

    input:
        path vcf                   // 定相后的 VCF 文件
        path vcf_tbi               // VCF 索引
        path alignment             // BAM 比对文件
        path alignment_index       // 比对文件索引
        path fasta                 // 参考基因组 FASTA
        path fasta_fai             // 参考基因组索引 (.fai)
        val sample_id              // 样本标识符
        val output_dir             // 输出目录

    output:
        path "${sample_id}.tagged.bam", emit: alignment
        path "${sample_id}.tagged.bam.bai", emit: index

    when:
    output_dir != 'NO_OUTPUT'

    script:
    """
    # 运行 WhatsHap 单倍型标记
    whatshap haplotag \\
        --reference ${fasta} \\
        --output "${sample_id}.tagged.bam" \\
        ${vcf} \\
        ${alignment}

    # 创建索引
    samtools index "${sample_id}.tagged.bam"
    """
}

// ============================================================================
// WHATSHAP_STATS - 定相统计
// ============================================================================
process WHATSHAP_STATS {
    tag "STATS on ${vcf.baseName}"
    label 'process_low'
    label 'whatshap'

    input:
        path vcf                   // 定相后的 VCF 文件
        path vcf_tbi               // VCF 索引 (可选，用于区域统计)
        val sample_id              // 样本标识符
        val output_dir             // 输出目录

    output:
        path "${sample_id}.phasing.stats.txt", emit: stats
        path "${sample_id}.phasing.stats.json", emit: json

    when:
    output_dir != 'NO_OUTPUT'

    script:
    """
    # 生成定相统计报告
    whatshap stats \\
        --output "${sample_id}.phasing.stats.txt" \\
        --json "${sample_id}.phasing.stats.json" \\
        ${vcf}
    """
}