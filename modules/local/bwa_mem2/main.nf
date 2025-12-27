/*
 * BWA 比对模块 (智能选择版本)
 * 
 * 功能：将 PE 测序数据比对到参考基因组，排序后输出 CRAM/BAM 格式并建立索引
 * 工具：bwa-mem2 (内存 >= 64GB) 或 bwa mem (内存 < 64GB) + samtools
 * 资源：由 modules.config 定义 (默认 16 CPU, 64GB 内存)
 * 
 * 智能决策：根据实际分配的内存自动选择比对工具
 *   - >= 64GB: 使用 bwa-mem2 (更快)
 *   - <  64GB: 使用 bwa mem (内存友好)
 * 
 * 索引要求：bwa/bwa-mem2 索引文件需与 fasta 同目录同前缀
 */
process BWA_MEM2 {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)  // meta: 样本元信息; reads: [R1.fq.gz, R2.fq.gz]
    path  fasta                   // 参考基因组 fasta (索引文件需同目录同前缀)
    val   output_format           // 输出格式: 'cram' (默认) 或 'bam'

    output:
    tuple val(meta), path("*.{cram,bam}"), path("*.{crai,bai}"), emit: alignment
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = meta.read_group ?: "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA"
    def format = output_format ?: 'cram'
    
    // 内存阈值: 64GB - 在 Groovy 层面判断并选择 aligner
    def mem_bytes = task.memory ? task.memory.toBytes() : 0
    def mem_threshold = 64L * 1024L * 1024L * 1024L
    def aligner = mem_bytes >= mem_threshold ? 'bwa-mem2' : 'bwa'
    
    def output_ext = format == 'cram' ? 'cram' : 'bam'
    def index_ext = format == 'cram' ? 'crai' : 'bai'
    def sort_opts = format == 'cram' ? "--reference ${fasta} -O cram" : "-O bam"
    
    """
    echo "使用比对工具: ${aligner} (分配内存: ${task.memory})"

    ${aligner} mem \\
        -t ${task.cpus} \\
        -R "${read_group}" \\
        ${args} \\
        ${fasta} \\
        ${reads[0]} \\
        ${reads[1]} \\
        | samtools sort -@ ${task.cpus} ${sort_opts} -o ${prefix}.${output_ext} -

    samtools index -@ ${task.cpus} ${prefix}.${output_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aligner: ${aligner}
        ${aligner}: \$(${aligner} version 2>&1 | head -n1)
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
