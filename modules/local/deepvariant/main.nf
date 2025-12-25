/*
 * DeepVariant 模块
 * 
 * 功能：基于深度学习的变异检测
 * 工具：Google DeepVariant
 * 
 * 说明：DeepVariant 使用 GPU 加速效果更好，但也支持 CPU 模式
 */
process DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(alignment), path(alignment_index)  // BAM/CRAM 及其索引 (.bai/.crai)
    path  fasta
    path  fasta_fai
    val   model_type    // 模型类型: 'WGS', 'WES' (默认), 'PACBIO', 'ONT_R104', 'HYBRID_PACBIO_ILLUMINA'
    path  intervals     // 可选，目标区域 BED 文件

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf      // 变异结果
    tuple val(meta), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), emit: gvcf, optional: true  // gVCF (可选)
    tuple val(meta), path("*.visual_report.html")          , emit: report   // 可视化报告
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def model = model_type ?: 'WES'  // 默认使用 WES 模型
    def regions_cmd = intervals ? "--regions=${intervals}" : ''
    def gvcf_cmd = params.deepvariant_gvcf ?: false ? "--output_gvcf=${prefix}.g.vcf.gz" : ''
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=${model} \\
        --ref=${fasta} \\
        --reads=${alignment} \\
        --output_vcf=${prefix}.vcf.gz \\
        ${gvcf_cmd} \\
        ${regions_cmd} \\
        --num_shards=${task.cpus} \\
        --intermediate_results_dir=./tmp \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(cat /opt/deepvariant/bin/VERSION || echo "unknown")
    END_VERSIONS
    """
}
