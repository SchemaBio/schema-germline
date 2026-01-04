/*
 * Stranger - STR Annotation Tool
 *
 * 功能：注释来自 ExpansionHunter 和 TRGT 的输出文件
 * 工具：stranger
 *
 * 说明：
 *   - 根据重复序列大小标注病理含义
 *   - 添加 STR_STATUS 字段 (normal/pre_mutation/full_mutation)
 *   - 支持 GRCh37/GRCh38 基因组版本
 */

process STRANGER {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // ExpansionHunter 输出 VCF
    path  repeats_file                          // 重复定义 JSON 文件
    val   family_id                             // 家族 ID (可选)
    val   is_trgt                               // 是否为 TRGT 输出

    output:
    tuple val(meta), path("*.stranger.vcf.gz"), path("*.stranger.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def repeats_cmd = repeats_file.name != 'NO_FILE' ? "--repeats-file ${repeats_file}" : ''
    def family_cmd = family_id ? "--family-id ${family_id}" : ''
    def trgt_cmd = is_trgt ? '--trgt' : ''
    """
    stranger \\
        ${vcf} \\
        ${repeats_cmd} \\
        ${family_cmd} \\
        ${trgt_cmd} \\
        ${args} | bgzip -c > ${prefix}.stranger.vcf.gz

    tabix -f -p vcf ${prefix}.stranger.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stranger: \$(stranger --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * STRANGER_MULTISAMPLE - 多样本 STR 注释
 *
 * 功能：同时注释多样本 VCF 文件
 */
process STRANGER_MULTISAMPLE {
    tag "stranger_multi"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 多样本 VCF
    path  repeats_file                          // 重复定义 JSON 文件

    output:
    tuple val(meta), path("*.stranger.vcf.gz"), path("*.stranger.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def repeats_cmd = repeats_file.name != 'NO_FILE' ? "--repeats-file ${repeats_file}" : ''
    """
    stranger \\
        ${vcf} \\
        ${repeats_cmd} \\
        ${args} | bgzip -c > ${prefix}.stranger.vcf.gz

    tabix -f -p vcf ${prefix}.stranger.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stranger: \$(stranger --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}
