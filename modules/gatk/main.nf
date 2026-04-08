// GATK 变异检测模块
// 用途：标记重复、合并比对、质控统计、变异检测
// 包含的process：
//   MARKDUPLICATES - 标记PCR重复
//   COLLECTQCMETRICS - 收集质控指标
//   LEFTALIGNANDTRIMVARIANTS - 左对齐和标准化变异
//   MUTECT2 - 线粒体分析

process MARKDUPLICATES {
    tag "MARKDUPLICATES on $sample_id"
    label 'gatk'
    label 'process_medium'
    publishDir "${params.output}/02.Alignment", mode: 'copy'

    input:
        path alignment      // BAM/CRAM 文件 (文件名包含 sample_id)
        path alignment_index // BAM/CRAM 索引文件
        val fasta           // 参考基因组路径

    output:
        path("*markdup.{cram,bam}"), emit: alignment
        path("*.markdup.metrics.txt"), emit: metrics

    script:
    // 从文件名提取 sample_id
    def sample_id = alignment.baseName.replaceAll(/\.(bam|cram)$/, '')
    """
    gatk MarkDuplicates \\
        -I ${alignment} \\
        -O ${sample_id}.markdup.${alignment.name.endsWith('.bam') ? 'bam' : 'cram'} \\
        -M ${sample_id}.markdup.metrics.txt \\
        --CREATE_INDEX false \\
        --REFERENCE_SEQUENCE ${fasta}
    """
}