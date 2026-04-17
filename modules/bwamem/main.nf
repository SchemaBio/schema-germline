// BWA MEM 比对模块
// 用途：将测序reads比对到参考基因组
// 用法：
//   - 输入：样本ID、FASTQ文件、参考基因组、Read Group ID
//   - 输出：BAM文件及索引
//   - 支持 BWA-MEM 和 BWA-MEM2 两种算法

process BWAMEM {
    tag "BWAMEM on $sample_id"
    label 'process_high'
    label 'mapping'

    input:
        tuple val(sample_id), path(reads)
        path fasta
        path fasta_index
        val rgid

    output:
        path("${sample_id}.bam"), emit: bam
        path("${sample_id}.bam.bai"), emit: bai

    script:
    def threads = task.cpus
    """
    bwa mem -t ${threads} -M -R "@RG\\tID:${rgid}\\tSM:${sample_id}\\tPL:SCHEMABIO\\tPU:Somatic" $fasta $reads \\
        | samtools sort -@ ${threads} -O bam -o ${sample_id}.bam -
    samtools index -@ ${threads} ${sample_id}.bam
    """
}

process BWAMEM2 {
    tag "BWAMEM2 on $sample_id"
    label 'process_high'
    label 'mapping'

    input:
        tuple val(sample_id), path(reads)
        path fasta
        path fasta_index
        val rgid

    output:
        path("${sample_id}.bam"), emit: bam
        path("${sample_id}.bam.bai"), emit: bai

    script:
    def threads = task.cpus
    """
    bwa-mem2 mem -t ${threads} -M -R "@RG\\tID:${rgid}\\tSM:${sample_id}\\tPL:SCHEMABIO\\tPU:Germline" $fasta $reads \\
        | samtools sort -@ ${threads} -O bam -o ${sample_id}.bam -
    samtools index -@ ${threads} ${sample_id}.bam
    """
}
