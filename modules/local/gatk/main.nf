/*
 * GATK 模块集合
 * 
 * 包含：CreateSequenceDictionary, MarkDuplicates, CollectQCMetrics, 
 *       LeftAlignAndTrimVariants, Mutect2 (线粒体模式)
 */

/*
 * CreateSequenceDictionary - 创建参考基因组字典文件
 */
process GATK_CREATESEQUENCEDICTIONARY {
    tag "$fasta"
    label 'process_low'

    input:
    path fasta

    output:
    path "*.dict"      , emit: dict
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = fasta.baseName
    """
    gatk CreateSequenceDictionary \\
        --REFERENCE ${fasta} \\
        --OUTPUT ${prefix}.dict \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP 'GATK v\\K[0-9.]+')
    END_VERSIONS
    """
}

/*
 * MarkDuplicates - 标记/去除 PCR 重复
 *
 * 注意：GATK 4.6.x 输出 CRAM 时不自动生成索引，改用 samtools index
 */
process GATK_MARKDUPLICATES {
    tag "$meta.id"

    input:
    tuple val(meta), path(alignment), path(alignment_index)  // BAM/CRAM 及其索引 (.bai/.crai)
    val   fasta_path      // CRAM 输出时需要参考序列路径
    val   output_format   // 输出格式: 'cram' (默认) 或 'bam'

    output:
    tuple val(meta), path("*.md.{cram,bam}"), emit: alignment
    tuple val(meta), path("*.metrics.txt")                                   , emit: metrics
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = output_format ?: 'cram'
    def output_file = format == 'cram' ? "${prefix}.md.cram" : "${prefix}.md.bam"
    def ref_cmd = format == 'cram' ? "--REFERENCE_SEQUENCE ${fasta_path}" : ''
    def index_cmd = format == 'cram' ? "samtools index -@ ${task.cpus} ${output_file}" : ''
    """
    gatk MarkDuplicates \\
        --INPUT ${alignment} \\
        --OUTPUT ${output_file} \\
        --METRICS_FILE ${prefix}.metrics.txt \\
        --CREATE_INDEX false \\
        ${ref_cmd} \\
        ${args}

    # 使用 samtools 生成 CRAM 索引
    if [ "${format}" == "cram" ]; then
        samtools index -@ ${task.cpus} ${output_file}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP 'GATK v\\K[0-9.]+')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

/*
 * CollectMultipleMetrics - 收集多种 QC 指标
 */
process GATK_COLLECTQCMETRICS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(alignment), path(alignment_index)  // BAM/CRAM 及其索引 (.bai/.crai)
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("*.metrics*"), emit: metrics  // 各类 QC 指标文件
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk CollectMultipleMetrics \\
        --INPUT ${alignment} \\
        --OUTPUT ${prefix}.metrics \\
        --REFERENCE_SEQUENCE ${fasta} \\
        --PROGRAM CollectAlignmentSummaryMetrics \\
        --PROGRAM CollectInsertSizeMetrics \\
        --PROGRAM QualityScoreDistribution \\
        --PROGRAM CollectGcBiasMetrics \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP 'GATK v\\K[0-9.]+')
    END_VERSIONS
    """
}

/*
 * LeftAlignAndTrimVariants - 左对齐并修剪变异
 */
process GATK_LEFTALIGNANDTRIMVARIANTS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)
    path  fasta
    path  fasta_fai
    path  dict

    output:
    tuple val(meta), path("*.normalized.vcf.gz"), path("*.normalized.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk LeftAlignAndTrimVariants \\
        --variant ${vcf} \\
        --reference ${fasta} \\
        --output ${prefix}.normalized.vcf.gz \\
        --split-multi-allelics \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP 'GATK v\\K[0-9.]+')
    END_VERSIONS
    """
}

/*
 * Mutect2 - 线粒体变异检测模式
 *
 * 针对线粒体特点优化：
 *   - 启用 --mitochondria-mode (自动调整参数适配 mtDNA)
 *   - 不使用 tumor/normal 配对模式
 *   - 支持低频异质性变异检测
 *
 * 注意：使用主参考基因组（如GRCh38），无需单独指定mt参考
 */
process GATK_MUTECT2_MT {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(alignment), path(alignment_index)  // BAM/CRAM 及其索引 (.bai/.crai)
    path  fasta          // 参考基因组（包含线粒体序列）
    path  fasta_fai      // 参考基因组索引

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.vcf.gz.stats")                , emit: stats
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk Mutect2 \\
        --reference ${fasta} \\
        --input ${alignment} \\
        --output ${prefix}.mt.vcf.gz \\
        --mitochondria-mode \\
        --max-reads-per-alignment-start 0 \\
        --max-mnp-distance 0 \\
        -L MT,chrM,chrMT,M \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP 'GATK v\\K[0-9.]+')
    END_VERSIONS
    """
}
