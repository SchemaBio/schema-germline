// DEEPVARIANT 模块
// 用途：基于深度学习的 SNV/Indel 变异检测 (全外显子测序)
// 工具：Google DeepVariant (CPU 模式)
// 输出：gVCF (用于 cohort 合并) + VCF (单样本调用)
//
// 说明：
//   - 本模块固定使用 WES 模型，适用于全外显子测序数据
//   - intervals 为必填参数，必须提供捕获区域 BED 文件
//   - num_shards 控制并行分片数，默认使用 task.cpus

process DEEPVARIANT {
    tag "DEEPVARIANT on ${alignment.baseName}"
    label 'process_high'
    label 'deepvariant'

    input:
        path alignment           // BAM 比对文件
        path alignment_index     // 比对文件索引 (.bai)
        path fasta               // 参考基因组 FASTA
        path fasta_fai           // 参考基因组索引 (.fai)
        path fasta_dict          // 参考基因组字典 (.dict)
        path intervals           // 目标区域 BED 文件 (必填)
        val num_shards           // 并行分片数 (默认使用 task.cpus)

    output:
        path "*.g.vcf.gz", emit: gvcf
        path "*.g.vcf.gz.tbi", emit: gvcf_tbi
        path "*.vcf.gz", emit: vcf
        path "*.vcf.gz.tbi", emit: vcf_tbi
        path "*.visual_report.html", emit: report, optional: true

    script:
    def shards = num_shards ?: task.cpus
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam)$/, '')
    """
    # 设置参考基因组环境变量 (DeepVariant 需要)
    export REFERENCE_GENOME=${fasta}

    # 处理 BED 文件：对每个区域前后拓展 20bp
    # 这有助于捕获目标区域边界附近的变异
    awk 'BEGIN {OFS="\\t"} {
        if (\$1 ~ /^#/) { print; next }  # 保留注释行
        start = \$2 - 20
        end = \$3 + 20
        if (start < 0) start = 0
        print \$1, start, end, \$4, \$5, \$6
    }' ${intervals} > intervals_extended.bed

    # 运行 DeepVariant (CPU 模式，WES 模型)
    # --num_shards: 并行处理分片数，加速处理
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type WES \\
        --ref ${fasta} \\
        --reads ${alignment} \\
        --output_vcf ${sample_id}.vcf.gz \\
        --output_gvcf ${sample_id}.g.vcf.gz \\
        --regions intervals_extended.bed \\
        --num_shards ${shards}

    # 创建索引 (DeepVariant 会自动创建，但确保存在)
    if [ ! -f ${sample_id}.vcf.gz.tbi ]; then
        tabix -p vcf ${sample_id}.vcf.gz
    fi
    if [ ! -f ${sample_id}.g.vcf.gz.tbi ]; then
        tabix -p vcf ${sample_id}.g.vcf.gz
    fi

    # 清理临时文件
    rm -f intervals_extended.bed
    """
}