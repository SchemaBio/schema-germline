/*
 * VEP (Variant Effect Predictor) 模块
 *
 * 功能：变异功能注释
 * 工具：Ensembl VEP
 *
 * 说明：
 *   - 支持 RefSeq 转录本注释
 *   - 支持 ClinVar、InterVar、cytoBand 等内置自定义注释
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
    path  clinvar_vcf                           // ClinVar VCF (可选)
    path  clinvar_tbi
    path  intervar_vcf                          // InterVar VCF (可选)
    path  intervar_tbi
    path  cytoband_bed                          // cytoBand BED (可选)
    path  cytoband_tbi
    path  dbnsfp_db                             // dbNSFP 数据库 (可选)
    path  dbnsfp_tbi
    path  spliceai_snv                          // SpliceAI SNV (可选)
    path  spliceai_snv_tbi
    path  spliceai_indel                        // SpliceAI Indel (可选)
    path  spliceai_indel_tbi
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

    // dbNSFP 字段
    def dbnsfp_fields = "rs_dbSNP,gnomAD_exomes_AF,gnomAD_exomes_EAS_AF,gnomAD_exomes_POPMAX_AF,gnomAD_exomes_AC,gnomAD_exomes_nhomalt,ExAC_AC,ExAC_EAS_AF,1000Gp3_AF,1000Gp3_EAS_AF,REVEL_score,REVEL_rankscore,M-CAP_score,M-CAP_rankscore,M-CAP_pred,GERP++_RS,GERP++_RS_rankscore,MVP_score,MVP_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP470way_mammalian,phyloP470way_mammalian_rankscore"

    // SpliceAI 字段
    def spliceai_fields = "SpliceAI_pred_DP_AG,SpliceAI_pred_DP_AL,SpliceAI_pred_DP_DG,SpliceAI_pred_DP_DL,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL"

    // 构建 custom 注释命令
    def clinvar_cmd = clinvar_vcf.name != 'NO_FILE' ? "--custom file=${clinvar_vcf},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNHGVS" : ''
    def intervar_cmd = intervar_vcf.name != 'NO_FILE' ? "--custom file=${intervar_vcf},short_name=InterVar,format=vcf,type=exact,coords=0,fields=Intervar" : ''
    def cytoband_cmd = cytoband_bed.name != 'NO_FILE' ? "--custom file=${cytoband_bed},short_name=cytoBand,format=bed,type=overlap,coords=0" : ''

    // 构建插件命令
    def dbnsfp_cmd = dbnsfp_db.name != 'NO_FILE' ? "--plugin dbNSFP,${dbnsfp_db},${dbnsfp_fields}" : ''
    def spliceai_cmd = (spliceai_snv.name != 'NO_FILE' && spliceai_indel.name != 'NO_FILE') ? "--plugin SpliceAI,snv=${spliceai_snv},indel=${spliceai_indel}" : ''

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
        ${cytoband_cmd} \\
        ${dbnsfp_cmd} \\
        ${spliceai_cmd} \\
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
