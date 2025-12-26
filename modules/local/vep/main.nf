/*
 * VEP (Variant Effect Predictor) 模块
 * 
 * 功能：变异功能注释
 * 工具：Ensembl VEP
 * 
 * 说明：支持多种注释数据库和插件，输出 VCF 或 JSON/TSV 格式
 */
process VEP {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 输入 VCF 及其索引
    path  fasta                                 // 参考基因组 (可选，用于 HGVS 注释)
    path  fasta_fai
    path  cache                                 // VEP 缓存目录
    path  plugins                               // VEP 插件目录 (可选)
    val   genome_assembly                       // 基因组版本: 'GRCh37' 或 'GRCh38'

    output:
    tuple val(meta), path("*.vep.vcf.gz"), path("*.vep.vcf.gz.tbi"), emit: vcf, optional: true
    tuple val(meta), path("*.vep.json.gz")                         , emit: json, optional: true
    tuple val(meta), path("*.vep.tsv")                             , emit: tsv, optional: true
    tuple val(meta), path("*.vep_summary.html")                    , emit: report
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assembly = genome_assembly ?: 'GRCh38'
    def fasta_cmd = fasta ? "--fasta ${fasta}" : ''
    def plugins_cmd = plugins ? "--dir_plugins ${plugins}" : ''
    def output_format = params.vep_output_format ?: 'vcf'  // vcf, json, tab
    def output_file = output_format == 'vcf' ? "${prefix}.vep.vcf.gz" : 
                      output_format == 'json' ? "${prefix}.vep.json.gz" : 
                      "${prefix}.vep.tsv"
    def compress_cmd = output_format in ['vcf', 'json'] ? '--compress_output bgzip' : ''
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${output_file} \\
        --${output_format} \\
        ${compress_cmd} \\
        --cache \\
        --dir_cache ${cache} \\
        --species homo_sapiens \\
        --assembly ${assembly} \\
        ${fasta_cmd} \\
        ${plugins_cmd} \\
        --offline \\
        --fork ${task.cpus} \\
        --stats_file ${prefix}.vep_summary.html \\
        --everything \\
        --force_overwrite \\
        ${args}

    # 为 VCF 输出创建索引
    if [[ "${output_format}" == "vcf" ]]; then
        tabix -p vcf ${output_file}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(vep --help 2>&1 | grep -oP 'ensembl-vep\\s+:\\s+\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}

/*
 * VEP_DOWNLOAD_CACHE - 下载 VEP 缓存数据
 * 
 * 用于首次运行或更新缓存
 */
process VEP_DOWNLOAD_CACHE {
    tag "$genome_assembly"
    label 'process_low'
    storeDir "${params.vep_cache_dir ?: 'vep_cache'}"

    input:
    val genome_assembly  // 'GRCh37' 或 'GRCh38'

    output:
    path "homo_sapiens*", emit: cache
    path "versions.yml" , emit: versions

    script:
    def assembly = genome_assembly ?: 'GRCh38'
    """
    vep_install \\
        --AUTO c \\
        --SPECIES homo_sapiens \\
        --ASSEMBLY ${assembly} \\
        --CACHEDIR . \\
        --NO_UPDATE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(vep --help 2>&1 | grep -oP 'ensembl-vep\\s+:\\s+\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}
