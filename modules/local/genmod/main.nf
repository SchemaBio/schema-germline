/*
 * GenMod - Genetic Model Annotation
 *
 * 功能：VCF 遗传模式注释工具
 * 工具：genmod
 *
 * 说明：
 *   - 注释遗传继承模式 (AD, AR, X-linked 等)
 *   - 复合杂合子检测
 *   - CADD 分数注释
 *   - 家系分析支持
 *   - 支持 VEP 注释结果
 */

// ============================================================================
// GENMOD_ANNOTATE - 基础注释
// ============================================================================
process GENMOD_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 输入 VCF
    path  ped_file                              // PED 家系文件 (可选)
    path  cadd_file                             // CADD 文件 (可选)
    val   annotate_cadd                         // 是否注释 CADD
    val   annotate_kg                           // 是否注释 1000G
    val   annotate_gnomad                       // 是否注释 gnomAD
    val   split_variants                        // 拆分多等位基因
    val   vep                                   // 支持 VEP 注释

    output:
    tuple val(meta), path("*.genmod_annot.vcf.gz"), path("*.genmod_annot.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ped_cmd = ped_file.name != 'NO_FILE' ? "--family_file ${ped_file}" : ''
    def cadd_cmd = cadd_file.name != 'NO_FILE' ? "--cadd_file ${cadd_file}" : ''
    def cadd_annot = annotate_cadd ? '--cadd' : ''
    def kg_annot = annotate_kg ? '--thousand_g' : ''
    def gnomad_annot = annotate_gnomad ? '--gnomad' : ''
    def split_cmd = split_variants ? '--split_variants' : ''
    def vep_cmd = vep ? '--vep' : ''
    """
    genmod annotate \
        ${vcf} \
        ${ped_cmd} \
        ${cadd_cmd} \
        ${cadd_annot} \
        ${kg_annot} \
        ${gnomad_annot} \
        ${split_cmd} \
        ${vep_cmd} \
        ${args} | bgzip -c > ${prefix}.genmod_annot.vcf.gz

    tabix -f -p vcf ${prefix}.genmod_annot.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmod: \$(genmod --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}

// ============================================================================
// GENMOD_MODELS - 遗传模式注释
// ============================================================================
process GENMOD_MODELS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 输入 VCF
    path  ped_file                              // PED 家系文件 (可选)
    val   split_variants                        // 拆分多等位基因
    val   phased                                // 输入已定相
    val   strict                                // 严格模式

    output:
    tuple val(meta), path("*.genmod_models.vcf.gz"), path("*.genmod_models.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ped_cmd = ped_file.name != 'NO_FILE' ? "--family_file ${ped_file}" : ''
    def split_cmd = split_variants ? '--split_variants' : ''
    def phased_cmd = phased ? '--phased' : ''
    def strict_cmd = strict ? '--strict' : ''
    """
    genmod models \
        ${vcf} \
        ${ped_cmd} \
        ${split_cmd} \
        ${phased_cmd} \
        ${strict_cmd} \
        ${args} | bgzip -c > ${prefix}.genmod_models.vcf.gz

    tabix -f -p vcf ${prefix}.genmod_models.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmod: \$(genmod --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}

// ============================================================================
// GENMOD_SORT - 按评分排序
// ============================================================================
process GENMOD_SORT {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 输入 VCF
    val   sort_by                              // 排序方式: 'score', 'position'

    output:
    tuple val(meta), path("*.genmod_sorted.vcf.gz"), path("*.genmod_sorted.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                                        , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sort_cmd = sort_by == 'score' ? '--score' : '--position'
    """
    genmod sort \
        ${vcf} \
        ${sort_cmd} \
        ${args} | bgzip -c > ${prefix}.genmod_sorted.vcf.gz

    tabix -f -p vcf ${prefix}.genmod_sorted.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmod: \$(genmod --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}

// ============================================================================
// GENMOD_SCORE - 变异评分
// ============================================================================
process GENMOD_SCORE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 输入 VCF
    path  ped_file                              // PED 家系文件 (可选)
    val   weight_cadd                           // CADD 权重
    val   weight_gnomad                         // gnomAD 权重
    val   weight_csq                            // CSQ 权重
    val   rank_threshold                        // 排名阈值

    output:
    tuple val(meta), path("*.genmod_score.vcf.gz"), path("*.genmod_score.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                                       , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ped_cmd = ped_file.name != 'NO_FILE' ? "--family_file ${ped_file}" : ''
    def cadd_w = weight_cadd != null ? "--cadd ${weight_cadd}" : ''
    def gnomad_w = weight_gnomad != null ? "--gnomad ${weight_gnomad}" : ''
    def csq_w = weight_csq != null ? "--csq ${weight_csq}" : ''
    def rank_cmd = rank_threshold != null ? "--rank_threshold ${rank_threshold}" : ''
    """
    genmod score \
        ${vcf} \
        ${ped_cmd} \
        ${cadd_w} \
        ${gnomad_w} \
        ${csq_w} \
        ${rank_cmd} \
        ${args} | bgzip -c > ${prefix}.genmod_score.vcf.gz

    tabix -f -p vcf ${prefix}.genmod_score.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmod: \$(genmod --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}

// ============================================================================
// GENMOD_FILTER - 变异过滤
// ============================================================================
process GENMOD_FILTER {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 输入 VCF
    val   filter_threshold                     // 过滤阈值
    val   remove_hom_germline                  // 去除同义突变
    val   filter_synonymous                    // 过滤同义突变
    val   minimum_score                        // 最低评分

    output:
    tuple val(meta), path("*.genmod_filter.vcf.gz"), path("*.genmod_filter.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                                        , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def thresh_cmd = filter_threshold != null ? "--filter_threshold ${filter_threshold}" : ''
    def syn_cmd = filter_synonymous ? '--filter_synonymous' : ''
    def hom_cmd = remove_hom_germline ? '--remove_hom_germline' : ''
    def min_cmd = minimum_score != null ? "--min_score ${minimum_score}" : ''
    """
    genmod filter \
        ${vcf} \
        ${thresh_cmd} \
        ${syn_cmd} \
        ${hom_cmd} \
        ${min_cmd} \
        ${args} | bgzip -c > ${prefix}.genmod_filter.vcf.gz

    tabix -f -p vcf ${prefix}.genmod_filter.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmod: \$(genmod --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}
