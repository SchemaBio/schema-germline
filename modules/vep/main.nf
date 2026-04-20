// VEP (Ensembl Variant Effect Predictor) 模块
// 用途：对 VCF 文件进行变异功能注释
// 工具：Ensembl VEP release 115.2
// 输出：注释后的 VCF 文件，包含基因、转录本、蛋白等注释信息
//
// 主要参数说明：
//   - cache: 使用本地缓存加速注释
//   - offline: 离线模式，不查询在线数据库
//   - pick: 选择每个变异最显著的注释结果
//   - everything: 输出所有可用注释字段
//   - fork: 多线程并行处理
//   - custom: 自定义数据库注释 (gnomAD, AlphaMissense, etc.)
//   - plugin: VEP 插件 (FlankingSequence, AnnotateClinVar, MissenseZscoreTranscript)
//
// ============================================================================
// VEP_ANNOTATE - SNV/Indel 变异注释
// ============================================================================
process VEP_ANNOTATE {
    tag "VEP on ${vcf.baseName}"
    label 'process_medium'
    label 'vep'

    input:
        path vcf               // 输入 VCF 文件
        path vcf_tbi           // VCF 索引文件 (.tbi)
        path fasta             // 参考基因组 FASTA
        path fasta_fai         // 参考基因组索引 (.fai)
        val genome_assembly    // 基因组版本: GRCh37 或 GRCh38
        val cache_dir          // VEP 缓存目录 (可选，默认使用内置缓存)
        val tag             // 输出文件前缀 (默认使用样本 ID)
        path vep_db_dir        // 自定义数据库根目录 (SchemaBio_Bundle)

    output:
        path "${sample_id}.${tag}.vep.vcf", emit: vep_vcf

    script:
    def assembly = genome_assembly ?: 'GRCh38'
    def db_prefix = assembly == 'GRCh37' ? 'hg19' : 'hg38'
    def sample_id = vcf.baseName.replaceAll(/\.(vcf|vcf\.gz)$/, '')
    def threads = Math.min(task.cpus as int, 8)

    // 设置 VEP 缓存目录
    def vep_cache = cache_dir ? cache_dir : '/opt/vep/.vep'
    def schema_bundle = '/schema_bundle'
    """
    cache_str="Uploaded_variation,Location,REF_ALLELE,Allele,Consequence,IMPACT,DOMAINS,Feature,DISTANCE,EXON,INTRON,SYMBOL,STRAND,HGNC_ID,HGVSc,HGVSp,HGVSg,MAX_AF,Protein_position,Amino_acids,Codons,PUBMED,Existing_variation"
    custom_str="cytoBand,CLNSIG,CLNDN,CLNSTAR,"
    schema_str=""
    self_plugin_str="FlankingSequence,MissenseZscore"

    vep \
        --offline --cache \
        --dir_cache ${cache_dir} --merged \
        --dir_plugins ${cache_dir}/Plugin \
        --force_overwrite --fork ${threads} \
        -i ${vcf} -o ${sample_id}.${tag}.vep.vcf \
        --format vcf --vcf \
        --fa ${fasta} \
        --shift_3prime 1 --assembly ${assembly} --no_escape --check_existing -exclude_predicted --uploaded_allele --show_ref_allele --numbers --domains \
        --total_length --hgvs --hgvsg --symbol --ccds --uniprot --max_af --pubmed \
        --transcript_filter "stable_id match N[MR]_" \
        --plugin AnnotateClinVar,clinvar_file=${schema_bundle}/${db_prefix}_clinvar.vcf.gz,fields=CLNSIG,CLNDN,CLNSTAR \
        --custom file=${schema_bundle}/${db_prefix}_cytoBand.bed.gz,short_name=cytoBand,format=bed,type=overlap,coords=0 \
        --custom file=${schema_bundle}/${db_prefix}_gnomad.v4.1.filtered.vcf.gz,short_name=GnomAD,format=vcf,type=exact,coords=0,fields=Exomes_AC_XY%Genomes_AC_XY%Gnomad_AC_XY \
        --custom file=${schema_bundle}/${db_prefix}_pangolin.vcf.gz,short_name=Pangolin,format=vcf,type=exact,coords=0,fields=Exomes_AC_XY%Genomes_AC_XY%Gnomad_AC_XY \
        --custom file=${schema_bundle}/${db_prefix}_EVOScore2.vcf.gz,short_name=EVOScore2,format=vcf,type=exact,coords=0,fields=Exomes_AC_XY%Genomes_AC_XY%Gnomad_AC_XY \
        --custom file=${schema_bundle}/${db_prefix}_AlphaMissense.v3.vcf.gz,short_name=AlphaMissense,format=vcf,type=exact,coords=0,fields=Exomes_AC_XY%Genomes_AC_XY%Gnomad_AC_XY \
        --plugin FlankingSequence,10 \
        --plugin MissenseZscoreTranscript,${schema_bundle}/missenseByTranscript.hg38.v4.1.bed \
        --fields "${cache_str},${custom_str},${schema_str},${self_plugin_str}"
    """
}

