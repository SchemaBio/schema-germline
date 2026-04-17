// Parabricks GPU加速工具模块
// 用途：使用NVIDIA Clara Parabricks进行GPU加速的生物信息学分析
// 包含的process：
//   PB_FQ2BAM - GPU加速的BWA MEM比对
//   PB_DEEPVARIANT - GPU加速的DeepVariant变异检测


process PB_FQ2BAM {
    tag "PB_FQ2BAM on $sample_id"
    label 'process_medium'
    label 'parabricks'

    input:
        tuple val(sample_id), path(reads)
        path fasta
        path fasta_index
        val rgid

    output:
        path("${sample_id}.marked.bam"), emit: bam
        path("${sample_id}.marked.bam.bai"), emit: bai

    script:
    def threads = task.cpus
    def gpu_devices = task.gpu ? task.gpu : 'all'
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
}

process PB_DEEPVARIANT {
    tag "PB_DEEPVARIANT on ${alignment.baseName}"
    label 'process_medium'
    label 'parabricks'

    input:
        path alignment           // BAM 比对文件
        path alignment_index     // 比对文件索引 (.bai)
        path fasta               // 参考基因组 FASTA
        path fasta_fai           // 参考基因组索引 (.fai)
        path fasta_dict          // 参考基因组字典 (.dict)
        path intervals           // 目标区域 BED 文件 (可选)
        val model_type           // 模型类型: WGS/WES/PACBIO/ONT_R104 (默认 WES)
        val num_shards           // 并行分片数 (默认使用 task.cpus)

    output:
        path "*.g.vcf.gz", emit: gvcf
        path "*.g.vcf.gz.tbi", emit: gvcf_tbi
        path "*.vcf.gz", emit: vcf
        path "*.vcf.gz.tbi", emit: vcf_tbi

    script:
    def model = model_type ?: 'WES'
    def shards = num_shards ?: task.cpus
    def gpu_devices = task.gpu ? task.gpu : 'all'
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam)$/, '')
    """
    # 处理 BED 文件：对每个区域前后拓展 20bp
    if [ "${intervals}" != "NO_FILE" ] && [ -s "${intervals}" ]; then
        awk 'BEGIN {OFS="\\t"} {
            if (\$1 ~ /^#/) { print; next }
            start = \$2 - 20
            end = \$3 + 20
            if (start < 0) start = 0
            print \$1, start, end
        }' ${intervals} > intervals_extended.bed
        REGIONS_PARAM="--intervals intervals_extended.bed"
    else
        REGIONS_PARAM=""
    fi

    # 运行 Parabricks DeepVariant (GPU 加速)
    # --model-type: 选择预训练模型
    # --gvcf: 输出 gVCF 格式
    pbrun deepvariant \\
        --ref ${fasta} \\
        --in-bam ${alignment} \\
        --out-vcf ${sample_id}.vcf \\
        --model-type ${model} \\
        --gvcf ${sample_id}.g.vcf \\
        \${REGIONS_PARAM} \\
        --num-workers ${shards} \\
        --gpu-devices ${gpu_devices}

    # 压缩并索引 VCF
    bgzip -c ${sample_id}.vcf > ${sample_id}.vcf.gz
    tabix -p vcf ${sample_id}.vcf.gz

    bgzip -c ${sample_id}.g.vcf > ${sample_id}.g.vcf.gz
    tabix -p vcf ${sample_id}.g.vcf.gz

    # 清理临时文件
    rm -f intervals_extended.bed ${sample_id}.vcf ${sample_id}.g.vcf
    """
}