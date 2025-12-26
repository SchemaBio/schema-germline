/*
 * ExpansionHunter 模块
 * 
 * 功能：短串联重复序列 (STR) 扩展检测
 * 工具：Illumina ExpansionHunter
 * 
 * 说明：用于检测与遗传病相关的 STR 扩展，如亨廷顿病、脆性X综合征等
 */
process EXPANSIONHUNTER {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(alignment), path(alignment_index)  // BAM/CRAM 及其索引
    path  fasta
    path  fasta_fai
    path  variant_catalog                                     // STR 位点目录 JSON 文件

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.json")                        , emit: json
    tuple val(meta), path("*_realigned.bam")               , emit: bam, optional: true
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sex = meta.sex ?: 'female'  // 默认 female，影响 X 染色体 STR 分析
    def realign_cmd = params.expansionhunter_output_bam ?: false ? "--output-prefix ${prefix}" : ''
    """
    ExpansionHunter \\
        --reads ${alignment} \\
        --reference ${fasta} \\
        --variant-catalog ${variant_catalog} \\
        --output-prefix ${prefix} \\
        --sex ${sex} \\
        --threads ${task.cpus} \\
        ${realign_cmd} \\
        ${args}

    # 压缩并索引 VCF
    bgzip -c ${prefix}.vcf > ${prefix}.vcf.gz
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$(ExpansionHunter --version 2>&1 | grep -oP 'ExpansionHunter v\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}

/*
 * EXPANSIONHUNTER_DENOVO - 无先验 STR 发现模式
 * 
 * 用于发现新的 STR 扩展位点，不依赖已知目录
 */
process EXPANSIONHUNTER_DENOVO {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(alignment), path(alignment_index)
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("*.str_profile.json"), emit: profile
    tuple val(meta), path("*.locus.tsv")       , emit: locus
    tuple val(meta), path("*.motif.tsv")       , emit: motif
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ExpansionHunterDenovo profile \\
        --reads ${alignment} \\
        --reference ${fasta} \\
        --output-prefix ${prefix} \\
        --min-anchor-mapq 50 \\
        --max-irr-mapq 40 \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter_denovo: \$(ExpansionHunterDenovo --version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * EXPANSIONHUNTER_DENOVO_MERGE - 合并多样本 STR profile
 * 
 * 用于群体分析，合并多个样本的 profile 进行比较
 */
process EXPANSIONHUNTER_DENOVO_MERGE {
    tag "merge"
    label 'process_low'

    input:
    path profiles        // 多个 .str_profile.json 文件
    path manifest        // 样本清单 TSV (sample_id, case/control, profile_path)
    path fasta
    path fasta_fai

    output:
    path "merged.multisample_profile.json", emit: merged_profile
    path "outlier_locus.tsv"              , emit: outlier_locus
    path "outlier_motif.tsv"              , emit: outlier_motif
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ExpansionHunterDenovo merge \\
        --reference ${fasta} \\
        --manifest ${manifest} \\
        --output-prefix merged \\
        ${args}

    # 异常值检测
    ExpansionHunterDenovo outlier \\
        --manifest ${manifest} \\
        --multisample-profile merged.multisample_profile.json \\
        --output outlier_locus.tsv \\
        --locus

    ExpansionHunterDenovo outlier \\
        --manifest ${manifest} \\
        --multisample-profile merged.multisample_profile.json \\
        --output outlier_motif.tsv \\
        --motif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter_denovo: \$(ExpansionHunterDenovo --version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}
