/*
 * BWA 比对模块 (智能选择版本)
 * 
 * 功能：将 PE 测序数据比对到参考基因组，排序后输出 CRAM/BAM 格式并建立索引
 * 工具：bwa-mem2 (内存 >= 32GB) 或 bwa mem (内存 < 32GB) + samtools
 * 
 * 智能决策：根据分配的内存自动选择比对工具
 *   - >= 32GB: 使用 bwa-mem2 (更快)
 *   - <  32GB: 使用 bwa mem (内存友好)
 */
process BWA_MEM2 {
    tag "$meta.id"       // 任务标签，显示样本ID便于追踪
    label 'process_high' // 资源标签，在 nextflow.config 中定义 CPU/内存

    input:
    tuple val(meta), path(reads)  // meta: 样本元信息 (id, read_group等); reads: [R1.fq.gz, R2.fq.gz]
    path  bwamem2_index            // bwa-mem2 索引目录
    path  bwa_index                // bwa 索引目录 (降级时使用)
    path  fasta                    // 参考基因组 fasta，用于 CRAM 压缩
    val   output_format            // 输出格式: 'cram' (默认) 或 'bam'

    output:
    tuple val(meta), path("*.{cram,bam}"), path("*.{crai,bai}"), emit: alignment  // 比对结果 + 索引
    path "versions.yml"                                        , emit: versions   // 软件版本信息

    when:
    task.ext.when == null || task.ext.when  // 条件执行控制

    script:
    def args = task.ext.args ?: ''  // 额外参数，可在 modules.config 中通过 ext.args 传入
    def prefix = task.ext.prefix ?: "${meta.id}"  // 输出文件前缀
    def read_group = meta.read_group ?: "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA"  // 默认 read group
    def format = output_format ?: 'cram'  // 默认输出 CRAM
    
    // 内存阈值: 32GB (单位: bytes)
    def mem_threshold = 32 * 1024 * 1024 * 1024
    def use_bwamem2 = task.memory?.toBytes() >= mem_threshold
    
    // 根据内存选择比对工具和索引
    def aligner = use_bwamem2 ? 'bwa-mem2' : 'bwa'
    def index_path = use_bwamem2 ? bwamem2_index : bwa_index
    def index_suffix = use_bwamem2 ? '.0123' : '.bwt'
    
    if (format == 'cram') {
        """
        # 自动查找索引前缀
        INDEX=\$(find -L ${index_path} -name "*${index_suffix}" | sed 's/${index_suffix}\$//')

        echo "使用比对工具: ${aligner} (分配内存: ${task.memory})"

        # 比对 -> 排序 -> 输出 CRAM
        ${aligner} mem \\
            -t ${task.cpus} \\
            -R "${read_group}" \\
            ${args} \\
            \$INDEX \\
            ${reads[0]} \\
            ${reads[1]} \\
            | samtools sort -@ ${task.cpus} --reference ${fasta} -O cram -o ${prefix}.cram -

        # 建立 CRAM 索引
        samtools index -@ ${task.cpus} ${prefix}.cram

        # 记录软件版本
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            aligner: ${aligner}
            ${aligner}: \$(${aligner} version 2>&1 | head -n1)
            samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
        END_VERSIONS
        """
    } else {
        """
        # 自动查找索引前缀
        INDEX=\$(find -L ${index_path} -name "*${index_suffix}" | sed 's/${index_suffix}\$//')

        echo "使用比对工具: ${aligner} (分配内存: ${task.memory})"

        # 比对 -> 排序 -> 输出 BAM
        ${aligner} mem \\
            -t ${task.cpus} \\
            -R "${read_group}" \\
            ${args} \\
            \$INDEX \\
            ${reads[0]} \\
            ${reads[1]} \\
            | samtools sort -@ ${task.cpus} -O bam -o ${prefix}.bam -

        # 建立 BAM 索引
        samtools index -@ ${task.cpus} ${prefix}.bam

        # 记录软件版本
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            aligner: ${aligner}
            ${aligner}: \$(${aligner} version 2>&1 | head -n1)
            samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
        END_VERSIONS
        """
    }
}
