// FASTP 质控模块
// 用途：对原始测序数据进行质量控制和接头去除
// 用法：
//   - 输入：样本ID、双端FASTQ文件
//   - 输出：过滤后的FASTQ、JSON和HTML质控报告
//   - 自动检测并去除接头序列

process FASTP {
    tag "FASTP on $sample_id"
    label 'process_medium'
    label 'mapping'
    publishDir "${params.output}/01.QC", mode: 'copy'

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id), path("${sample_id}.clean_*.fq.gz"), emit: clean_reads
        path("${sample_id}.fastp.stat.json"), emit: json_report
        path("${sample_id}.fastp.stat.html"), emit: html_report

    script:
    def threads = Math.min(task.cpus as int, 16)  // fastp 最多支持 16 线程
    """
    fastp \\
        -i ${read1} \\
        -I ${read2} \\
        -o ${sample_id}.clean_1.fq.gz \\
        -O ${sample_id}.clean_2.fq.gz \\
        -w ${threads} \\
        -j ${sample_id}.fastp.stat.json \\
        -h ${sample_id}.fastp.stat.html \\
        --detect_adapter_for_pe
    """
}
