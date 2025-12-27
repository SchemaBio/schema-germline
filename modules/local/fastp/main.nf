/*
 * FASTP 质控过滤模块
 * 
 * 功能：对 PE 测序数据进行质量控制和过滤
 * 工具：fastp
 * 
 * 输出：
 *   - 过滤后的 fastq 文件
 *   - JSON/HTML 质控报告
 */
process FASTP {
    tag "$meta.id"
    label 'process_low'  // fastp 内存需求低

    input:
    tuple val(meta), path(reads)  // meta: 样本元信息; reads: [R1.fq.gz, R2.fq.gz]

    output:
    tuple val(meta), path("*_R{1,2}.fq.gz"), emit: reads    // 过滤后的 fastq
    tuple val(meta), path("*.json")        , emit: json     // JSON 报告
    tuple val(meta), path("*.html")        , emit: html     // HTML 报告
    path "versions.yml"                    , emit: versions // 软件版本信息

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''  // 额外参数，可在 modules.config 中配置
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = Math.min(task.cpus as int, 16)  // fastp 最多支持 16 线程
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${prefix}_R1.fq.gz \\
        --out2 ${prefix}_R2.fq.gz \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --thread ${threads} \\
        --detect_adapter_for_pe \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """
}
