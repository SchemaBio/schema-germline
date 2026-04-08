// Samtools 模块
// 用途：为 BAM/CRAM 文件创建索引
// 用法：
//   - 输入：BAM/CRAM 文件
//   - 输出：索引文件 (.bai/.crai)

process SAMTOOLS_INDEX {
    tag "SAMTOOLS_INDEX on ${alignment.baseName}"
    label 'process_low'
    publishDir "${params.output}/02.Alignment", mode: 'copy'

    input:
        path alignment

    output:
        path alignment, emit: alignment
        path "*.{bai,crai}", emit: index

    script:
    def threads = Math.min(task.cpus as int, 4)
    """
    samtools index -@ ${threads} ${alignment}
    """
}