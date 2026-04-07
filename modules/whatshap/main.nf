/*
 * WhatsHap 单倍型 phasing 模块
 *
 * 功能：基于比对 reads 信息对 VCF 进行单倍型分相 (phasing)
 * 工具：WhatsHap
 *
 * 说明：
 *   - 使用 reads 比对信息推断等位基因的顺式关系
 *   - 支持 indels phasing
 *   - 输出带 HP 标签的 phased VCF
 */

process WHATSHAP_PHASE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 输入 VCF 及其索引
    tuple val(meta2), path(bam), path(bai)     // BAM/CRAM 及其索引
    path  fasta                                 // 参考基因组
    path  fasta_fai

    output:
    tuple val(meta), path("*.phase.vcf.gz"), path("*.phase.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_cmd = bam.name.endsWith('.cram') ? "--reference ${fasta}" : ''
    """
    # 运行 WhatsHap phasing
    whatshap phase \\
        --indels \\
        --reference ${fasta} \\
        --output ${prefix}.phase.vcf \\
        ${args} \\
        ${vcf} \\
        ${bam}

    # 压缩并索引
    bgzip -f ${prefix}.phase.vcf
    tabix -f -p vcf ${prefix}.phase.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version 2>&1 | sed 's/whatshap //')
    END_VERSIONS
    """
}

/*
 * WHATSHAP_GENOTYPE - 基因型质量提升
 *
 * 功能：使用 reads 信息提升基因型质量分数 (GT, GQ)
 * 工具：WhatsHap
 */
process WHATSHAP_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)
    tuple val(meta2), path(bam), path(bai)
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("*.genotyped.vcf.gz"), path("*.genotyped.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    whatshap genotype \\
        --reference ${fasta} \\
        --output ${prefix}.genotyped.vcf \\
        ${args} \\
        ${vcf} \\
        ${bam}

    bgzip -f ${prefix}.genotyped.vcf
    tabix -f -p vcf ${prefix}.genotyped.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version 2>&1 | sed 's/whatshap //')
    END_VERSIONS
    """
}

/*
 * WHATSHAP_STATS - Phasing 统计
 *
 * 功能：生成 phasing 结果的统计报告
 */
process WHATSHAP_STATS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(phase_vcf), path(phase_vcf_tbi)
    tuple val(meta2), path(bam), path(bai)

    output:
    tuple val(meta), path("*.stats.txt"), emit: stats
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    whatshap stats \\
        ${phase_vcf} \\
        > ${prefix}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version 2>&1 | sed 's/whatshap //')
    END_VERSIONS
    """
}
