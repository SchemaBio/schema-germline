// BWA MEM 比对模块
// 用途：将测序reads比对到参考基因组
// 用法：
//   - 输入：样本ID、FASTQ文件、参考基因组、输出格式(bam/cram)、Read Group ID
//   - 输出：BAM/CRAM文件及索引
//   - 支持 BWA-MEM 和 BWA-MEM2 两种算法

process BWAMEM {
    tag "BWAMEM on $sample_id"
    label 'process_high'
    label 'mapping'

    input:
        tuple val(sample_id), path(reads)
        path fasta
        path fasta_index
        val output_format
        val rgid

    output:
        path("${sample_id}.cram"), emit: cram, optional: true
        path("${sample_id}.cram.crai"), emit: crai, optional: true
        path("${sample_id}.bam"), emit: bam, optional: true
        path("${sample_id}.bam.bai"), emit: bai, optional: true

    script:
    def threads = task.cpus
    if (output_format == 'bam') {
        """
        bwa mem -t ${threads} -M -R "@RG\\tID:${rgid}\\tSM:${sample_id}\\tPL:SCHEMABIO\\tPU:Somatic" $fasta $reads \\
            | samtools sort -@ ${threads} -O bam -o ${sample_id}.bam -
        samtools index -@ ${threads} ${sample_id}.bam
        """
    } else {
        """
        bwa mem -t ${threads} -M -R "@RG\\tID:${rgid}\\tSM:${sample_id}\\tPL:SCHEMABIO\\tPU:Somatic" $fasta $reads \\
            | samtools sort -@ ${threads} \\
            --reference ${fasta} -O cram --output-fmt-option version=3.0 -o ${sample_id}.cram -
        samtools index -@ ${threads} ${sample_id}.cram
        """
    }
}

process BWAMEM2 {
    tag "BWAMEM2 on $sample_id"
    label 'process_high'
    label 'mapping'

    input:
        tuple val(sample_id), path(reads)
        path fasta
        path fasta_index
        val output_format
        val rgid

    output:
        path("${sample_id}.cram"), emit: cram, optional: true
        path("${sample_id}.cram.crai"), emit: crai, optional: true
        path("${sample_id}.bam"), emit: bam, optional: true
        path("${sample_id}.bam.bai"), emit: bai, optional: true

    script:
    def threads = task.cpus
    if (output_format == 'bam') {
        """
        bwa-mem2 mem -t ${threads} -M -R "@RG\\tID:${rgid}\\tSM:${sample_id}\\tPL:SCHEMABIO\\tPU:Somatic" $fasta $reads \\
            | samtools sort -@ ${threads} -O bam -o ${sample_id}.bam -
        samtools index -@ ${threads} ${sample_id}.bam
        """
    } else {
        """
        bwa-mem2 mem -t ${threads} -M -R "@RG\\tID:${rgid}\\tSM:${sample_id}\\tPL:SCHEMABIO\\tPU:Somatic" $fasta $reads \\
            | samtools sort -@ ${threads} \\
            --reference ${fasta} -O cram --output-fmt-option version=3.0 -o ${sample_id}.cram -
        samtools index -@ ${threads} ${sample_id}.cram
        """
    }
}
