// BWA MEM 比对模块
// 用途：将测序 reads 比对到参考基因组
// 用法：
//   - 输入：样本元信息、FASTQ 文件、参考基因组、BWA 索引文件集合
//   - 输出：BAM 文件
//   - 支持 BWA-MEM 和 BWA-MEM2 两种算法
//
// 输入格式：
//   - meta: [id: sample_id, read_group: "@RG\\tID:...\\tSM:..."]
//   - reads: [read1.fastq.gz, read2.fastq.gz]
//   - bwa_indices: *.amb, *.ann, *.bwt, *.pac, *.sa 文件集合

process BWAMEM {
    tag "BWAMEM on ${meta.id}"
    label 'process_high'
    label 'mapping'

    input:
        tuple val(meta), path(reads)  // [meta, [read1, read2]]
        path fasta
        path bwa_indices              // BWA 索引文件集合 (*.amb, *.ann, *.bwt, *.pac, *.sa)

    output:
        tuple val(meta), path("${meta.id}.bam"), emit: bam

    script:
    def sample_id = meta.id
    def rg = meta.read_group ?: "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:SchemaBio:\\tPU:Germline"
    def threads = task.cpus
    """
    bwa mem \\
        -t ${threads} \\
        -M \\
        -R "${rg}" \\
        ${fasta} \\
        ${reads[0]} ${reads[1]} \\
        | samtools sort -@ ${threads} -O bam -o ${sample_id}.bam -
    samtools index -@ ${threads} ${sample_id}.bam
    """
}

process BWAMEM2 {
    tag "BWAMEM2 on ${meta.id}"
    label 'process_huge'
    label 'mapping'

    input:
        tuple val(meta), path(reads)  // [meta, [read1, read2]]
        path fasta
        path bwamem2_indices          // BWA-MEM2 索引文件集合 (*.0123, *.bwt.2bit.64)

    output:
        tuple val(meta), path("${meta.id}.bam"), emit: bam

    script:
    def sample_id = meta.id
    def rg = meta.read_group ?: "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:SchemaBio:\\tPU:Germline"
    def threads = task.cpus
    """
    bwa-mem2 mem \\
        -t ${threads} \\
        -M \\
        -R "${rg}" \\
        ${fasta} \\
        ${reads[0]} ${reads[1]} \\
        | samtools sort -@ ${threads} -O bam - -o ${sample_id}.bam
    samtools index -@ ${threads} ${sample_id}.bam
    """
}
