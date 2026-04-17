// EXPANSIONHUNTER 模块
// 用途：短 tandem repeat (STR) 扩展检测
// 工具：ExpansionHunter
// 输出：STR 扩展检测结果 VCF + JSON
//
// 说明：
//   - 专门用于检测三核苷酸重复扩展疾病 (如亨廷顿病、脆性X综合征等)
//   - 需要提供 STR 位点定义文件 (variant catalog, JSON 格式)
//   - 支持性别特异性分析 (影响 X/Y 染色体 STR 的检测)
//   - variant catalog 可从 ExpansionHunter 官方获取，或自定义

// ============================================================================
// EXPANSIONHUNTER - STR 扩展检测
// ============================================================================
process EXPANSIONHUNTER {
    tag "EXPANSIONHUNTER on ${alignment.baseName}"
    label 'process_low'
    label 'expansionhunter'

    input:
        path alignment           // BAM 比对文件
        path alignment_index     // 比对文件索引
        path fasta               // 参考基因组 FASTA
        path fasta_fai           // 参考基因组索引 (.fai)
        path variant_catalog     // STR 位点定义文件 (JSON 格式)
        val sex                  // 样本性别: 'male', 'female', 或 'unknown' (默认 unknown)
        val min_anchor           // 最小锚定长度 (默认 8)
        val max_irr_mapping      // 最大不完美匹配距离 (默认 100)
        val output_dir           // 输出目录

    output:
        path "${sample_id}.vcf.gz", emit: vcf
        path "${sample_id}.vcf.gz.tbi", emit: vcf_tbi
        path "${sample_id}.json", emit: json
        path "${sample_id}_reads.json", emit: reads_json, optional: true

    when:
    output_dir != 'NO_OUTPUT'

    script:
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam)$/, '')
    def sex_param = sex ? "--sex ${sex}" : ''
    def anchor = min_anchor ?: 8
    def irr = max_irr_mapping ?: 100
    """
    # 运行 ExpansionHunter
    expansionhunter \\
        --reads ${alignment} \\
        --reference ${fasta} \\
        --variant-catalog ${variant_catalog} \\
        --output-prefix ${sample_id} \\
        --min-anchor ${anchor} \\
        --max-irr-mapping ${irr} \\
        ${sex_param} \\
        --output-dir .

    # 压缩 VCF 并创建索引
    bgzip -c "${sample_id}.vcf" > "${sample_id}.vcf.gz"
    tabix -p vcf "${sample_id}.vcf.gz"

    # 清理原始 VCF
    rm -f "${sample_id}.vcf"
    """
}