// Samtools 模块
// 用途：为 BAM/CRAM 文件创建索引、性别检测
// 包含的process：
//   SAMTOOLS_INDEX - 创建索引
//   SEX_CHECK_SRY - 基于 SRY 读数的性别检测

process SAMTOOLS_INDEX {
    tag "SAMTOOLS_INDEX on ${alignment.baseName}"
    label 'process_low'
    label 'mapping'
    publishDir "${params.output}/02.Alignment", mode: 'copy'

    input:
        path alignment

    output:
        path alignment, emit: alignment
        path "*.{bai,crai}", emit: index

    script:
    def threads = Math.min(task.cpus as int, 4)
    """
    samtools index -@ ${threads} ${alignment}
    """
}

process SEX_CHECK_SRY {
    tag "SEX_CHECK_SRY on ${alignment.baseName}"
    label 'process_low'
    label 'mapping'
    publishDir "${params.output}/01.QC", mode: 'copy'

    input:
        path alignment
        path alignment_index
        val genome_assembly    // 'GRCh37' 或 'GRCh38'
        val threshold          // 判定男性的 SRY reads 阈值

    output:
        path "*.sex_check.json", emit: json

    script:
    // SRY 基因坐标 (包含上下游 500bp 缓冲区)
    // GRCh37: Y:2654896-2655792 -> Y:2654396-2656292
    // GRCh38: Y:2786855-2787741 -> Y:2786355-2788241
    def sry_regions = [
        'GRCh37': ['nochr': 'Y:2654396-2656292', 'chr': 'chrY:2654396-2656292'],
        'GRCh38': ['nochr': 'Y:2786355-2788241', 'chr': 'chrY:2786355-2788241']
    ]
    def assembly = genome_assembly ?: 'GRCh38'
    def sry_nochr = sry_regions[assembly]['nochr']
    def sry_chr = sry_regions[assembly]['chr']
    def thresh = threshold ?: 10
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam|cram)$/, '')
    """
    # 检测染色体命名方式
    CHR_STYLE=\$(samtools view -H ${alignment} 2>/dev/null | grep -E "^@SQ.*SN:(chr)?Y" | head -1 | grep -o "SN:chr" || echo "")

    if [ -n "\${CHR_STYLE}" ]; then
        SRY_REGION="${sry_chr}"
        CHR_PREFIX="chr"
    else
        SRY_REGION="${sry_nochr}"
        CHR_PREFIX="nochr"
    fi

    # 统计 SRY 区域的 reads 数 (排除未比对 reads 和重复)
    # -F 2052: 排除未比对的 reads (0x4) 和重复 reads (0x400)
    SRY_READS=\$(samtools view -c -F 2052 ${alignment} \${SRY_REGION} 2>/dev/null || echo "0")

    # 性别判定
    if [ "\${SRY_READS}" -ge ${thresh} ]; then
        INFERRED_SEX="male"
    else
        INFERRED_SEX="female"
    fi

    # 输出 JSON 格式结果
    cat <<EOF > ${sample_id}.sex_check.json
    {
        "sample_id": "${sample_id}",
        "genome_assembly": "${assembly}",
        "chr_style": "\${CHR_PREFIX}",
        "sry_region": "\${SRY_REGION}",
        "sry_reads": \${SRY_READS},
        "threshold": ${thresh},
        "inferred_sex": "\${INFERRED_SEX}"
    }
    EOF
    """
}