/*
 * UPD (Uniparental Disomy) 检测模块
 * 
 * 功能：检测单亲二倍体区域
 * 工具：UPD (https://github.com/bjhall/upd)
 * 
 * 说明：需要 trio 样本（proband + 父母）的 VCF 文件
 */
process UPD {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 多样本 VCF (trio)
    val   proband_id                            // 先证者样本 ID
    val   father_id                             // 父亲样本 ID
    val   mother_id                             // 母亲样本 ID

    output:
    tuple val(meta), path("*.upd.bed")        , emit: regions
    tuple val(meta), path("*.upd.sites.bed")  , emit: sites
    tuple val(meta), path("*.upd_summary.txt"), emit: summary, optional: true
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    upd \\
        --vcf ${vcf} \\
        --proband ${proband_id} \\
        --father ${father_id} \\
        --mother ${mother_id} \\
        --out-regions ${prefix}.upd.bed \\
        --out-sites ${prefix}.upd.sites.bed \\
        ${args}

    # 生成摘要统计
    if [ -s ${prefix}.upd.bed ]; then
        echo "UPD Summary for ${proband_id}" > ${prefix}.upd_summary.txt
        echo "================================" >> ${prefix}.upd_summary.txt
        echo "Total UPD regions: \$(wc -l < ${prefix}.upd.bed)" >> ${prefix}.upd_summary.txt
        echo "" >> ${prefix}.upd_summary.txt
        cat ${prefix}.upd.bed >> ${prefix}.upd_summary.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        upd: \$(upd --version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * UPD_SITES - 仅输出信息位点（不做区域合并）
 * 
 * 用于下游可视化或自定义分析
 */
process UPD_SITES {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)
    val   proband_id
    val   father_id
    val   mother_id

    output:
    tuple val(meta), path("*.informative_sites.bed"), emit: sites
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    upd \\
        --vcf ${vcf} \\
        --proband ${proband_id} \\
        --father ${father_id} \\
        --mother ${mother_id} \\
        --sites-only \\
        --out-sites ${prefix}.informative_sites.bed \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        upd: \$(upd --version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}
