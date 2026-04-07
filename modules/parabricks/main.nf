/*
 * ParaBricks 模块集合
 *
 * 功能：NVIDIA GPU 加速的基因组分析工具
 * 官方文档：https://docs.nvidia.com/clara/parabricks/
 *
 * 包含：
 *   - FQ2BAM: GPU 加速 BWA 比对 + 排序 + 标记重复 (替代 BWA + Samtools + GATK MarkDuplicates)
 *   - DEEPVARIANT: GPU 加速 DeepVariant 变异检测
 *   - HAPLOTYPECALLER: GPU 加速 GATK HaplotypeCaller
 *   - MUTECT2: GPU 加速 GATK Mutect2 (体细胞变异)
 *
 * 注意：
 *   - 需要 NVIDIA GPU 和 CUDA 环境
 *   - 需要指定 GPU 数量 (--gpu-devices 参数)
 *   - 输出格式与 CPU 版本兼容
 */

/*
 * FQ2BAM - GPU 加速 BWA 比对 + 排序 + 标记重复
 *
 * 功能：一步完成 BWA-MEM 比对、排序、标记重复
 * 输入：FASTQ 文件对 + 参考基因组 + BWA 索引
 * 输出：去重后的 BAM/CRAM 文件 + 索引 + 质控指标
 *
 * GPU 需求：建议 16GB+ 显存
 */
process PARABRICKS_FQ2BAM {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(reads)               // FASTQ 文件对
    path  fasta                                // 参考基因组 FASTA
    path  "${fasta}.fai"                       // FASTA 索引
    path  "${fasta}.dict"                      // FASTA 字典
    path  bwa_index                            // BWA 索引文件集合
    val   output_format                        // 输出格式: 'cram' (默认) 或 'bam'
    val   gpu_devices                          // GPU 设备 ID (如: 0 或 0,1,2,3)

    output:
    tuple val(meta), path("*.{cram,bam}")                , emit: alignment   // 比对结果
    tuple val(meta), path("*.{crai,bai}")               , emit: index       // 索引文件
    tuple val(meta), path("*.metrics.txt")              , emit: metrics     // 标记重复指标
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = meta.read_group ?: "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA"
    def format = output_format ?: 'cram'
    def output_ext = format == 'cram' ? 'cram' : 'bam'
    def gpus = gpu_devices ?: '0'

    """
    pbrun fq2bam \\
        --ref ${fasta} \\
        --in-fq ${reads[0]} ${reads[1]} \\
        --out-file ${prefix}.${output_ext} \\
        --out-metrics ${prefix}.metrics.txt \\
        --out-qc-metrics \\
        --RG "${read_group}" \\
        --gpu-devices ${gpus} \\
        --num-gpus 1 \\
        ${format == 'cram' ? '--output-format CRAM' : '--output-format BAM'} \\
        ${args}

    # ParaBricks 4.0+ 不自动生成索引，需要手动索引
    samtools index -@ ${task.cpus} ${prefix}.${output_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(pbrun --version 2>&1 | head -n1 || echo "unknown")
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

/*
 * DEEPVARIANT - GPU 加速 DeepVariant 变异检测
 *
 * 功能：使用 GPU 加速 DeepVariant SNP/INDEL 检测
 * 输入：BAM/CRAM + 参考基因组
 * 输出：VCF 变异结果
 *
 * GPU 需求：建议 16GB+ 显存，支持多 GPU 并行
 */
process PARABRICKS_DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(alignment), path(alignment_index)  // BAM/CRAM 及索引
    path  fasta                                              // 参考基因组
    path  "${fasta}.fai"                                     // FASTA 索引
    val   model_type                                         // 模型类型: 'WGS', 'WES'
    path  intervals                                          // 可选，目标区域 BED
    val   gpu_devices                                        // GPU 设备 ID

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def model = model_type ?: 'WES'
    def gpus = gpu_devices ?: '0'
    def regions_param = intervals ? "--regions ${intervals}" : ''

    """
    pbrun deepvariant \\
        --ref ${fasta} \\
        --in-bam ${alignment} \\
        --out-vcf ${prefix}.vcf.gz \\
        --model-type ${model} \\
        --gpu-devices ${gpus} \\
        ${regions_param} \\
        ${args}

    # 索引 VCF 文件
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(pbrun --version 2>&1 | head -n1 || echo "unknown")
        deepvariant: \$(cat /opt/deepvariant/bin/VERSION || echo "unknown")
    END_VERSIONS
    """
}

/*
 * HAPLOTYPECALLER - GPU 加速 GATK HaplotypeCaller
 *
 * 功能：使用 GPU 加速 GATK HaplotypeCaller 变异检测
 * 输入：BAM/CRAM + 参考基因组 + 目标区域 (可选)
 * 输出：GVCF 或 VCF
 *
 * GPU 需求：建议 16GB+ 显存
 */
process PARABRICKS_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(alignment), path(alignment_index)  // BAM/CRAM 及索引
    path  fasta                                              // 参考基因组
    path  "${fasta}.fai"                                     // FASTA 索引
    path  "${fasta}.dict"                                    // FASTA 字典
    path  intervals                                          // 可选，目标区域 BED
    val   emit_gvcf                                          // 是否输出 GVCF: true/false
    val   gpu_devices                                        // GPU 设备 ID

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), emit: gvcf, optional: true
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gpus = gpu_devices ?: '0'
    def regions_param = intervals ? "--intervals ${intervals}" : ''
    def gvcf_param = emit_gvcf ? "--emit-gvcf --gvcf ${prefix}.g.vcf.gz" : ''

    """
    pbrun haplotypecaller \\
        --ref ${fasta} \\
        --in-bam ${alignment} \\
        --out-vcf ${prefix}.vcf.gz \\
        ${regions_param} \\
        ${gvcf_param} \\
        --gpu-devices ${gpus} \\
        ${args}

    # 索引 VCF 文件
    tabix -p vcf ${prefix}.vcf.gz
    if [ -f "${prefix}.g.vcf.gz" ]; then
        tabix -p vcf ${prefix}.g.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(pbrun --version 2>&1 | head -n1 || echo "unknown")
        gatk: \$(gatk --version 2>&1 | grep -oP 'GATK v\\K[0-9.]+')
    END_VERSIONS
    """
}

/*
 * MUTECT2 - GPU 加速 GATK Mutect2 (体细胞变异检测)
 *
 * 功能：使用 GPU 加速 Mutect2 体细胞变异检测
 * 输入：Tumor BAM + Normal BAM (可选) + 参考基因组
 * 输出：体细胞变异 VCF
 *
 * GPU 需求：建议 16GB+ 显存
 */
process PARABRICKS_MUTECT2 {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai)        // Tumor BAM/CRAM
    tuple val(normal_meta), path(normal_bam), path(normal_bai), optional: true  // Normal BAM (可选)
    path  fasta                                              // 参考基因组
    path  "${fasta}.fai"                                     // FASTA 索引
    path  "${fasta}.dict"                                    // FASTA 字典
    path  intervals                                          // 可选，目标区域 BED
    val   gpu_devices                                        // GPU 设备 ID

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gpus = gpu_devices ?: '0'
    def regions_param = intervals ? "--intervals ${intervals}" : ''
    def normal_param = normal_bam ? "--normal ${normal_bam}" : ''

    """
    pbrun mutect2 \\
        --ref ${fasta} \\
        --in-tumor ${tumor_bam} \\
        ${normal_param} \\
        --out-vcf ${prefix}.vcf.gz \\
        ${regions_param} \\
        --gpu-devices ${gpus} \\
        ${args}

    # 索引 VCF 文件
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(pbrun --version 2>&1 | head -n1 || echo "unknown")
        gatk: \$(gatk --version 2>&1 | grep -oP 'GATK v\\K[0-9.]+')
    END_VERSIONS
    """
}

/*
 * GERMLINE - GPU 加速胚系变异检测 (完整流程)
 *
 * 功能：一步完成比对 + 标记重复 + HaplotypeCaller
 * 输入：FASTQ + 参考基因组 + BWA 索引
 * 输出：VCF 变异结果
 *
 * 适用于：胚系单样本 WGS/WES 分析
 */
process PARABRICKS_GERMLINE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(reads)               // FASTQ 文件对
    path  fasta                                // 参考基因组
    path  "${fasta}.fai"                       // FASTA 索引
    path  "${fasta}.dict"                      // FASTA 字典
    path  bwa_index                            // BWA 索引
    path  intervals                            // 可选，目标区域 BED
    val   gpu_devices                          // GPU 设备 ID

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = meta.read_group ?: "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA"
    def gpus = gpu_devices ?: '0'
    def regions_param = intervals ? "--intervals ${intervals}" : ''

    """
    pbrun germline \\
        --ref ${fasta} \\
        --in-fq ${reads[0]} ${reads[1]} \\
        --out-vcf ${prefix}.vcf.gz \\
        --RG "${read_group}" \\
        ${regions_param} \\
        --gpu-devices ${gpus} \\
        ${args}

    # 索引 VCF 文件
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(pbrun --version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}