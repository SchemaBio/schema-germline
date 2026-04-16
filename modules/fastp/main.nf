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


process COLLECTQCMETRICS {
    tag "COLLECTQCMETRICS on $sample_id"
    label 'gatk'
    label 'process_medium'

    input:
        tuple val(sample_id), path(cram), path(crai)
        val   fasta
        val   fasta_dict
        val   bed

    output:
        path("${sample_id}.hs_metrics.txt"), emit: hs_metrics
        path("${sample_id}.insert_size_metrics.txt"), emit: insert_size_metrics
        path("${sample_id}.insert_size_histogram.pdf"), emit: insert_size_histogram
        path("${sample_id}.est_lib_complex_metrics.txt"), emit: lib_metrics

    script:
    """
    gatk BedToIntervalList \\
        --VALIDATION_STRINGENCY SILENT \\
        --TMP_DIR tmp \\
        --SEQUENCE_DICTIONARY ${fasta_dict} \\
        --INPUT ${bed} \\
        --OUTPUT Bait.interval_list

    gatk CollectHsMetrics \\
        --VALIDATION_STRINGENCY SILENT \\
        --TMP_DIR tmp \\
        --BAIT_INTERVALS Bait.interval_list \\
        --TARGET_INTERVALS Bait.interval_list \\
        --INPUT ${cram} \\
        --OUTPUT ${sample_id}.hs_metrics.txt

    gatk CollectInsertSizeMetrics \\
        --VALIDATION_STRINGENCY SILENT \\
        --TMP_DIR tmp \\
        --INPUT ${cram} \\
        --OUTPUT ${sample_id}.insert_size_metrics.txt \\
        --Histogram_FILE ${sample_id}.insert_size_histogram.pdf

    gatk EstimateLibraryComplexity \\
        --MAX_RECORDS_IN_RAM 303942330 \\
        --VALIDATION_STRINGENCY SILENT \\
        --TMP_DIR tmp \\
        --INPUT ${cram} \\
        --OUTPUT ${sample_id}.est_lib_complex_metrics.txt
    """
}
