/*
 * Samtools 模块集合
 *
 * 包含：Index
 */

/*
 * SAMTOOLS_INDEX - 为 BAM/CRAM 文件创建索引
 */
process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(alignment)  // BAM 或 CRAM 文件

    output:
    tuple val(meta), path(alignment), path("*.{bai,crai}"), emit: alignment
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools index \\
        -@ ${task.cpus} \\
        ${args} \\
        ${alignment}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
