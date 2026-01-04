/*
 * Slivar - Genetic Variant Filtering & Annotation
 *
 * 功能：使用 JavaScript 表达式进行 VCF 变异过滤和注释
 * 工具：slivar
 *
 * 说明：
 *   - 支持复杂过滤表达式
 *   - 支持 gnomAD 等数据库快速注释
 *   - 支持家系 (trio/family) 分析
 *   - 支持 de novo、复合杂合子检测
 */

process SLIVAR_EXPR {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 输入 VCF 及其索引
    path  ped_file                              // PED 系谱文件 (可选)
    path  gnotate                               // gnomad gnotate 文件 (可选)
    path  slivar_js                             // slivar-functions.js (可选)
    val   info_expr                             // INFO 字段过滤表达式
    val   sample_expr                           // 样本表达式
    val   trio_expr                             // 三联体表达式
    val   family_expr                           // 家族表达式
    val   pass_only                             // 只输出通过的变体

    output:
    tuple val(meta), path("*.slivar.vcf.gz"), path("*.slivar.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ped_cmd = ped_file.name != 'NO_FILE' ? "--ped ${ped_file}" : ''
    def gnotate_cmd = gnotate.name != 'NO_FILE' ? "-g ${gnotate}" : ''
    def js_cmd = slivar_js.name != 'NO_FILE' ? "--js ${slivar_js}" : ''
    def info_cmd = info_expr ? "--info '${info_expr}'" : ''
    def sample_cmd = sample_expr ? "--sample-expr '${sample_expr}'" : ''
    def trio_cmd = trio_expr ? "--trio '${trio_expr}'" : ''
    def family_cmd = family_expr ? "--family-expr '${family_expr}'" : ''
    def pass_cmd = pass_only ? '--pass-only' : ''
    """
    slivar expr \\
        --vcf ${vcf} \\
        ${ped_cmd} \\
        ${gnotate_cmd} \\
        ${js_cmd} \\
        ${info_cmd} \\
        ${sample_cmd} \\
        ${trio_cmd} \\
        ${family_cmd} \\
        ${pass_cmd} \\
        --out-vcf ${prefix}.slivar.vcf.gz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: \$(slivar version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * SLIVAR_COMPOUND_HETS - 复合杂合子检测
 *
 * 功能：检测复合杂合子变异（两个不同位点的杂合变异）
 */
process SLIVAR_COMPOUND_HETS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)
    path  ped_file                              // PED 文件 (可选)
    path  gnotate                               // gnomad gnotate 文件 (可选)

    output:
    tuple val(meta), path("*.compound.vcf.gz"), path("*.compound.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ped_cmd = ped_file.name != 'NO_FILE' ? "--ped ${ped_file}" : ''
    def gnotate_cmd = gnotate.name != 'NO_FILE' ? "-g ${gnotate}" : ''
    """
    slivar compound-hets \\
        --vcf ${vcf} \\
        ${ped_cmd} \\
        ${gnotate_cmd} \\
        --out-vcf ${prefix}.compound.vcf.gz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: \$(slivar version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * SLIVAR_TSV - 导出 TSV 格式
 *
 * 功能：将 VCF 导出为 TSV 表格便于分析
 */
process SLIVAR_TSV {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)
    path  ped_file                              // PED 文件 (可选)

    output:
    tuple val(meta), path("*.slivar.tsv.gz"), emit: tsv
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ped_cmd = ped_file.name != 'NO_FILE' ? "--ped ${ped_file}" : ''
    """
    slivar tsv \\
        --vcf ${vcf} \\
        ${ped_cmd} \\
        --output ${prefix}.slivar.tsv.gz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: \$(slivar version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * SLIVAR_DUO_DEL - 亲子缺失检测
 *
 * 功能：检测样本中的亲子缺失变异
 */
process SLIVAR_DUO_DEL {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)
    path  ped_file

    output:
    tuple val(meta), path("*.duo_del.tsv.gz"), emit: tsv
    path "versions.yml"                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    slivar duo-del \\
        --vcf ${vcf} \\
        --ped ${ped_file} \\
        --output ${prefix}.duo_del.tsv.gz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: \$(slivar version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}
