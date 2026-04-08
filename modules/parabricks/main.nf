// Parabricks GPU加速工具模块
// 用途：使用NVIDIA Clara Parabricks进行GPU加速的生物信息学分析
// 用法：
//   PB_BAM2FQ - 将BAM/CRAM文件转换为FASTQ格式
//   PB_DEEPVARIANT - GPU加速的DeepVariant变异检测
//   PB_FQ2BAM - GPU加速的BWA MEM比对


process PB_FQ2BAM {
    tag "PB_FQ2BAM on $sample_id"
    label 'process_medium'
    label 'parabricks'
    publishDir "${params.output}/02.Alignment", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)
        path fasta
        path fasta_index
        val output_format
        val rgid

    output:
        path("${sample_id}.marked.cram"), emit: cram, optional: true
        path("${sample_id}.marked.cram.crai"), emit: crai, optional: true
        path("${sample_id}.marked.bam"), emit: bam, optional: true
        path("${sample_id}.marked.bam.bai"), emit: bai, optional: true

    script:
    def threads = task.cpus
    def gpu_devices = task.gpu ? task.gpu : 'all'
    if (output_format == 'bam') {
        """
        pbrun fq2bam \\
            --ref ${fasta} \\
            --in-fq ${reads} \\
            --out-bam ${sample_id}.bam \\
            --rg-id "${rgid}" \\
            --rg-sm "${sample_id}" \\
            --rg-pl "SCHEMABIO" \\
            --rg-pu "Germline" \\
            --threads ${threads} \\
            --gpu-devices ${gpu_devices} \\
            --gpusort

        samtools index -@ ${threads} ${sample_id}.bam
        """
    } else {
        """
        pbrun fq2bam \\
            --ref ${fasta} \\
            --in-fq ${reads} \\
            --out-bam ${sample_id}.cram \\
            --rg-id "${rgid}" \\
            --rg-sm "${sample_id}" \\
            --rg-pl "SCHEMABIO" \\
            --rg-pu "Germline" \\
            --threads ${threads} \\
            --gpu-devices ${gpu_devices} \\
            --gpusort

        samtools index -@ ${threads} ${sample_id}.cram
        """
    }
}