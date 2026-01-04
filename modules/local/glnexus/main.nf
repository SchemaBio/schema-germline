/*
 * GLnexus 模块
 *
 * 功能：合并多个 gVCF 文件进行联合基因分型
 * 工具：GLnexus
 *
 * 说明：GLnexus 是 DeepVariant 推荐的 gVCF 合并工具
 */

process GLNEXUS {
    tag "$meta.id"
    label 'process_high'

    input:
    val   meta              // 家系/批次元信息
    path  gvcfs             // 多个 gVCF 文件
    path  gvcf_tbis         // gVCF 索引文件
    val   config            // 配置: 'DeepVariantWES', 'DeepVariantWGS', 'DeepVariant_unfiltered'

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cfg = config ?: 'DeepVariantWES'
    """
    glnexus_cli \\
        --config ${cfg} \\
        --threads ${task.cpus} \\
        ${args} \\
        ${gvcfs} \\
        | bcftools view - \\
        | bgzip -c -@ ${task.cpus} \\
        > ${prefix}.joint.vcf.gz

    tabix -p vcf ${prefix}.joint.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus: \$(glnexus_cli --help 2>&1 | head -n1 | grep -oP 'v[0-9.]+' || echo "unknown")
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
