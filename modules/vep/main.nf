/*
 * VEP (Variant Effect Predictor) 模块
 *
 * 功能：变异功能注释
 * 工具：Ensembl VEP
 *
 * 说明：
 *   - 支持 RefSeq 转录本注释
 *   - 支持 ClinVar、InterVar、gnomAD、dbSNP 等自定义注释
 *   - 支持 AlphaMissense、EVOScore2、Pangolin 等插件
 *   - 支持 dbNSFP、SpliceAI 插件
 *   - 支持用户自定义 --custom 和 --plugin 参数
 *   - 支持 pick 模式 (只选择一个转录本) 或全部转录本
 */
process VEP {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 输入 VCF 及其索引
    path  fasta                                 // 参考基因组
    path  fasta_fai
    path  cache                                 // VEP 缓存目录
    path  plugins                               // VEP 插件目录
    // 自定义数据库 (VCF 格式)
    path  clinvar_vcf                           // ClinVar VCF (可选)
    path  clinvar_tbi
    path  intervar_vcf                          // InterVar VCF (可选)
    path  intervar_tbi
    path  gnomad_vcf                            // gnomAD VCF (可选)
    path  gnomad_tbi
    path  dbsnp_vcf                             // dbSNP VCF (可选)
    path  dbsnp_tbi
    path  alphamissense_db                      // AlphaMissense DB (可选)
    path  alphamissense_tbi
    path  evoscore2_db                          // EVOScore2 DB (可选)
    path  evoscore2_tbi
    path  pangolin_db                           // Pangolin DB (可选)
    path  pangolin_tbi
    path  cytoband_bed                          // cytoBand BED (可选)
    path  cytoband_tbi
    // 插件数据库
    path  extra_files                           // 用户自定义数据库文件 (可选，支持多个)
    val   extra_custom                          // 用户自定义 --custom 参数字符串
    val   extra_plugins                         // 用户自定义 --plugin 参数字符串
    val   genome_assembly                       // 基因组版本: 'GRCh37' 或 'GRCh38'

    output:
    tuple val(meta), path("*.vep.vcf.gz"), path("*.vep.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.vep_summary.html")                    , emit: report
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assembly = genome_assembly ?: 'GRCh38'

    // 是否使用 pick 模式 (只选择一个转录本)
    def use_pick = params.vep_pick ?: true
    def pick_cmd = use_pick ? '--pick --pick_order mane_plus_clinical,mane_select,refseq,canonical,ensembl --pick_allele' : ''

    // 是否只保留 RefSeq 转录本 (用于 merged 缓存)
    def refseq_only = params.vep_refseq_only ?: false
    def transcript_filter = refseq_only ? '--transcript_filter "stable_id match N[MR]_"' : ''

    // AlphaMissense 字段
    def alphamissense_fields = "AlphaMissense_score,AlphaMissense_pred"

    // EVOScore2 字段
    def evoscore2_fields = "EVOScore2_score,EVOScore2_pred"

    // Pangolin 字段
    def pangolin_fields = "Pangolin_score,Pangolin_pred,Pangolin_rankscore"

    // 构建 custom 注释命令 (VCF 格式数据库)
    def clinvar_cmd = clinvar_vcf.name != 'NO_FILE' ? "--custom file=${clinvar_vcf},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNHGVS" : ''
    def intervar_cmd = intervar_vcf.name != 'NO_FILE' ? "--custom file=${intervar_vcf},short_name=InterVar,format=vcf,type=exact,coords=0,fields=Intervar" : ''
    def gnomad_cmd = gnomad_vcf.name != 'NO_FILE' ? "--custom file=${gnomad_vcf},short_name=gnomAD,format=vcf,type=exact,coords=0,fields=AF%AF_eas%af_popmax" : ''
    def dbsnp_cmd = dbsnp_vcf.name != 'NO_FILE' ? "--custom file=${dbsnp_vcf},short_name=dbSNP,format=vcf,type=exact,coords=0,fields=RS%CLNID%CLNREVSTAT" : ''
    def cytoband_cmd = cytoband_bed.name != 'NO_FILE' ? "--custom file=${cytoband_bed},short_name=cytoBand,format=bed,type=overlap,coords=0" : ''

    // 构建插件命令
    def alphamissense_cmd = alphamissense_db.name != 'NO_FILE' ? "--plugin AlphaMissense,${alphamissense_db},${alphamissense_fields}" : ''
    def evoscore2_cmd = evoscore2_db.name != 'NO_FILE' ? "--plugin EVOScore2,${evoscore2_db},${evoscore2_fields}" : ''
    def pangolin_cmd = pangolin_db.name != 'NO_FILE' ? "--plugin Pangolin,${pangolin_db},${pangolin_fields}" : ''

    // 用户自定义参数
    def user_custom_cmd = extra_custom ?: ''
    def user_plugin_cmd = extra_plugins ?: ''
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${prefix}.vep.vcf.gz \\
        --format vcf \\
        --vcf \\
        --compress_output bgzip \\
        --offline \\
        --cache \\
        --dir_cache ${cache} \\
        --dir_plugins ${plugins} \\
        --refseq \\
        --fasta ${fasta} \\
        --assembly ${assembly} \\
        --force_overwrite \\
        --fork ${task.cpus} \\
        --stats_file ${prefix}.vep_summary.html \\
        --shift_3prime 1 \\
        --no_escape \\
        --check_existing \\
        --exclude_predicted \\
        --uploaded_allele \\
        --show_ref_allele \\
        --numbers \\
        --domains \\
        --total_length \\
        --hgvs \\
        --hgvsg \\
        --symbol \\
        --ccds \\
        --uniprot \\
        --max_af \\
        --pubmed \\
        ${pick_cmd} \\
        ${transcript_filter} \\
        ${clinvar_cmd} \\
        ${intervar_cmd} \\
        ${gnomad_cmd} \\
        ${dbsnp_cmd} \\
        ${cytoband_cmd} \\
        ${alphamissense_cmd} \\
        ${evoscore2_cmd} \\
        ${pangolin_cmd} \\
        ${user_custom_cmd} \\
        ${user_plugin_cmd} \\
        ${args}

    tabix -p vcf ${prefix}.vep.vcf.gz

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
