/*
 * BCFtools 模块集合
 *
 * 包含：ROH (Runs of Homozygosity) 检测
 */

/*
 * BCFTOOLS_ROH - ROH 检测
 *
 * 功能：检测基因组中的纯合子区域（Runs of Homozygosity）
 * 用途：近亲繁殖分析、家系研究、遗传病研究
 * 工具：bcftools roh
 */
process BCFTOOLS_ROH {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // VCF 及其索引
    path  fasta          // 参考基因组
    path  fasta_fai

    output:
    tuple val(meta), path("*.roh"), emit: roh           // ROH 结果文件
    tuple val(meta), path("*.roh_af0.15.hom") , emit: hom_samples  // 高ROH样本
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # 使用 GATK 风格的 AF 计算模式
    bcftools roh \
        --threads ${task.cpus} \
        --AF-dflt 0.001 \
        --output ${prefix}.roh \
        ${args} \
        ${vcf}

    # 提取高 ROH 区域样本 (AF < 0.15 的纯合变异)
    bcftools view -i 'GT="hom"' ${vcf} | \
        bcftools +setGT - -- -t q -n . -u > ${prefix}.roh_af0.15.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}

/*
 * BCFTOOLS_STATS - VCF 统计
 *
 * 功能：生成 VCF 文件的统计信息
 */
process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("*.vcfstats.txt"), emit: stats
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools stats \
        --threads ${task.cpus} \
        --reference ${fasta} \
        ${vcf} \
        > ${prefix}.vcfstats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
