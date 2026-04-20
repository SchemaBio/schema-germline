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
        val output_dir             // 输出目录

    output:
        path "${sample_id}.phased.vcf.gz", emit: vcf
        path "${sample_id}.phased.vcf.gz.tbi", emit: vcf_tbi
        path "${sample_id}.phasing.json", emit: json, optional: true

    when:
    output_dir != 'NO_OUTPUT'

    script:
    """
    # 运行 WhatsHap 定相
    whatshap phase \\
        --reference ${fasta} \\
        -o "${sample_id}.phased.vcf" \\
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
