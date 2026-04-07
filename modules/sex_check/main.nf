/*
 * 性别质控模块
 *
 * 功能：通过 SRY 基因区域的 reads 覆盖判断样本性别
 * 工具：samtools
 *
 * 说明：SRY 基因位于 Y 染色体，男性样本应有显著覆盖，女性样本几乎无覆盖
 *       支持 GRCh37 和 GRCh38 两种参考基因组
 */

// SRY 基因坐标 (包含上下游 500bp 缓冲区)
// GRCh37: Y:2654896-2655792
// GRCh38: Y:2786855-2787741
def SRY_REGIONS = [
    'GRCh37': [nochr: 'Y:2654396-2656292', chr: 'chrY:2654396-2656292'],
    'GRCh38': [nochr: 'Y:2786355-2788241', chr: 'chrY:2786355-2788241']
]

process SEX_CHECK {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(alignment), path(alignment_index)  // BAM/CRAM 及其索引
    val   genome_assembly                                     // 'GRCh37' 或 'GRCh38'

    output:
    tuple val(meta), path("*.sex_check.txt"), emit: result
    tuple val(meta), env(INFERRED_SEX)      , emit: sex
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assembly = genome_assembly ?: 'GRCh38'
    def regions = SRY_REGIONS[assembly] ?: SRY_REGIONS['GRCh38']
    def sry_region_nochr = regions.nochr
    def sry_region_chr = regions.chr
    def declared_sex = meta.sex ?: 'unknown'
    // 阈值：reads 数 >= 10 判定为男性
    def threshold = params.sex_check_threshold ?: 10
    """
    # 检测染色体命名方式
    CHR_STYLE=\$(samtools view -H ${alignment} 2>/dev/null | grep -E "^@SQ.*SN:(chr)?Y" | head -1 | grep -o "SN:chr" || echo "")

    if [ -n "\${CHR_STYLE}" ]; then
        SRY_REGION="${sry_region_chr}"
    else
        SRY_REGION="${sry_region_nochr}"
    fi

    # 统计 SRY 区域的 reads 数（排除未比对 reads 和重复）
    # -F 2052: 排除未比对的 reads (0x4) 和重复 reads (0x400)
    SRY_READS=\$(samtools view -c -F 2052 ${alignment} \${SRY_REGION} 2>/dev/null || echo "0")

    # 性别推断
    if [ "\${SRY_READS}" -ge ${threshold} ]; then
        INFERRED_SEX="male"
    else
        INFERRED_SEX="female"
    fi

    # 一致性检查
    DECLARED="${declared_sex}"
    if [ "\${DECLARED}" != "unknown" ] && [ "\${DECLARED}" != "\${INFERRED_SEX}" ]; then
        STATUS="MISMATCH"
    else
        STATUS="PASS"
    fi

    # 输出结果
    cat <<-EOF > ${prefix}.sex_check.txt
    sample_id: ${meta.id}
    genome_assembly: ${assembly}
    chr_style: \$([ -n "\${CHR_STYLE}" ] && echo "with_chr" || echo "no_chr")
    sry_region: \${SRY_REGION}
    sry_reads: \${SRY_READS}
    threshold: ${threshold}
    inferred_sex: \${INFERRED_SEX}
    declared_sex: ${declared_sex}
    status: \${STATUS}
    EOF

    # 导出环境变量供 Nextflow 捕获
    export INFERRED_SEX

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

/*
 * SEX_CHECK_MULTISAMPLE - 批量性别检查汇总
 */
process SEX_CHECK_SUMMARY {
    tag "summary"
    label 'process_low'

    input:
    path sex_check_files  // 多个 .sex_check.txt 文件

    output:
    path "sex_check_summary.tsv", emit: summary
    path "versions.yml"         , emit: versions

    script:
    """
    # 生成汇总表头
    echo -e "sample_id\\tsry_reads\\tinferred_sex\\tdeclared_sex\\tstatus" > sex_check_summary.tsv

    # 解析每个文件
    for f in ${sex_check_files}; do
        SAMPLE=\$(grep "sample_id:" \$f | cut -d' ' -f2)
        READS=\$(grep "sry_reads:" \$f | cut -d' ' -f2)
        INFERRED=\$(grep "inferred_sex:" \$f | cut -d' ' -f2)
        DECLARED=\$(grep "declared_sex:" \$f | cut -d' ' -f2)
        STATUS=\$(grep "status:" \$f | cut -d' ' -f2)
        echo -e "\${SAMPLE}\\t\${READS}\\t\${INFERRED}\\t\${DECLARED}\\t\${STATUS}" >> sex_check_summary.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1)
    END_VERSIONS
    """
}
