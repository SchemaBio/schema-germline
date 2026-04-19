/*
 * Schema Germline Pipeline - 主入口文件
 *
 * 整合三个 workflow：
 *   - WES_SINGLE: 单人全外显子测序分析
 *   - WES_TRIO_PIPELINE: Trio 家系全外显子测序分析
 *   - CNV_BASELINE: CNVkit 参考基线构建
 *
 * 使用方法：
 *   nextflow run main.nf --config <config.json> -profile <profile>
 *   nextflow run main.nf -entry WES_SINGLE --config examples/single.json -profile docker
 *   nextflow run main.nf -entry WES_TRIO --config examples/trio.json -profile docker
 *   nextflow run main.nf -entry CNV_BASELINE --config examples/baseline.json -profile docker
 */

nextflow.enable.dsl = 2

// ============================================================================
// Include Workflows
// ============================================================================
include { WES_SINGLE } from './workflows/single/main'
include { WES_TRIO_PIPELINE as WES_TRIO } from './workflows/trio/main'
include { CNV_BASELINE } from './workflows/baseline/main'

// ============================================================================
// 参数定义
// ============================================================================
params {
    config = null      // JSON 配置文件路径
    entry  = null      // 显式指定入口: WES_SINGLE, WES_TRIO, CNV_BASELINE
}

// ============================================================================
// 主 workflow - 根据配置自动分发
// ============================================================================
workflow {

    // 检查配置文件
    if (!params.config) {
        exit 1, "必须提供配置文件: --config <path/to/config.json>"
    }

    // 解析 JSON 配置
    def config_path = file(params.config)
    if (!config_path.exists()) {
        exit 1, "配置文件不存在: ${params.config}"
    }

    def config = new groovy.json.JsonSlurper().parseText(config_path.text)

    // 确定执行入口
    def entry_name = params.entry ?: detect_workflow_type(config)

    log.info "=========================================="
    log.info "Schema Germline Pipeline"
    log.info "=========================================="
    log.info "配置文件: ${params.config}"
    log.info "执行入口: ${entry_name}"
    log.info "输出目录: ${config.outdir ?: params.outdir}"
    log.info "=========================================="

    // 根据入口分发执行
    switch (entry_name) {
        case 'WES_SINGLE':
            run_wes_single(config)
            break
        case 'WES_TRIO':
            run_wes_trio(config)
            break
        case 'CNV_BASELINE':
            run_cnv_baseline(config)
            break
        default:
            exit 1, "未知的入口类型: ${entry_name}"
    }
}

// ============================================================================
// Workflow 类型检测
// ============================================================================
def detect_workflow_type(config) {
    // Trio: 有 family_id 和 samples 数组（含 role）
    if (config.family_id && config.samples) {
        def has_roles = config.samples.any { it.role }
        if (has_roles) {
            return 'WES_TRIO'
        }
    }

    // Single: 有单个样本的 sample_id 和 read1/read2
    if (config.sample_id && config.read1 && config.read2) {
        return 'WES_SINGLE'
    }

    // Baseline: 有 samples 数组（无 role）和 cnvkit 配置
    if (config.samples && !config.samples.any { it.role }) {
        return 'CNV_BASELINE'
    }

    exit 1, "无法从配置文件识别 workflow 类型，请使用 -entry 参数显式指定"
}

// ============================================================================
// WES_SINGLE 执行器
// ============================================================================
def run_wes_single(config) {

    // 参考基因组通道
    def fasta = config.reference.fasta
    def fasta_fai = "${fasta}.fai"
    def fasta_dict = fasta.replaceFirst(/\.fa(sta)?$/, '.dict')

    ch_fasta = Channel.fromPath(fasta, checkIfExists: true)
    ch_fasta_fai = Channel.fromPath(fasta_fai, checkIfExists: true)
    ch_fasta_dict = Channel.fromPath(fasta_dict, checkIfExists: true)

    // BED 文件通道
    ch_target_bed = Channel.fromPath(config.target_bed, checkIfExists: true)
    ch_mt_bed = Channel.fromPath(config.mt_bed, checkIfExists: true)
    ch_snp_positions = Channel.fromPath(config.snp_positions, checkIfExists: true)

    // CNVkit 文件通道
    ch_cnvkit_reference = Channel.fromPath(config.cnvkit_reference, checkIfExists: true)
    ch_cnvkit_access_bed = config.cnvkit_access_bed ?
        Channel.fromPath(config.cnvkit_access_bed, checkIfExists: true) : Channel.empty()

    // ExpansionHunter catalog 通道
    ch_expansionhunter_catalog = Channel.fromPath(config.expansionhunter_catalog, checkIfExists: true)

    // 样本 reads 通道
    ch_reads = Channel.value([
        config.sample_id,
        [file(config.read1), file(config.read2)]
    ])

    // Read Group ID
    def rgid = "@RG\\tID:${config.sample_id}\\tSM:${config.sample_id}\\tPL:ILLUMINA"

    // 输出目录
    def outdir = config.outdir ?: params.outdir

    // 执行 workflow
    WES_SINGLE(
        ch_reads,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        ch_target_bed,
        ch_mt_bed,
        ch_snp_positions,
        rgid,
        config.use_gpu ?: false,
        config.genome_assembly ?: 'GRCh38',
        config.sex_check_threshold ?: 10,
        config.baf_max_depth ?: 200,
        config.deepvariant_model ?: 'WES',
        config.deepvariant_shards,
        // CNVkit
        ch_cnvkit_reference,
        ch_cnvkit_access_bed,
        config.cnvkit_annotate ?: true,
        config.cnvkit_split_size ?: 5000,
        config.cnvkit_min_target_size ?: 10000,
        config.cnvkit_method ?: 'cbs',
        config.cnvkit_threshold,
        config.cnvkit_ploidy ?: 2,
        config.cnvkit_drop_outliers ?: true,
        "${outdir}/cnv",
        // ExpansionHunter
        ch_expansionhunter_catalog,
        config.expansionhunter_min_anchor ?: 8,
        config.expansionhunter_max_irr ?: 100,
        "${outdir}/str",
        // TIEA-WES
        config.tiea_min_support ?: 10,
        config.tiea_min_softclip ?: 36,
        config.tiea_min_mapq ?: 20,
        config.tiea_cluster_window ?: 10,
        config.tiea_threads ?: 4,
        // WhatsHap
        config.whatshap_chromosomes,
        config.whatshap_ignore_rg ?: true,
        config.whatshap_ref_conf ?: 20,
        "${outdir}/variants/snv_indel/phased",
        // VEP
        config.vep_use_pick ?: true,
        config.vep_use_refseq_only ?: false,
        config.vep_cache_dir,
        config.vep_extra_args,
        config.vep_mei_upstream_distance ?: 5000,
        config.vep_mei_downstream_distance ?: 5000,
        config.vep_db_dir,
        config.vep_use_gnomad ?: true,
        config.vep_use_alphamissense ?: true,
        config.vep_use_evoscore ?: true,
        config.vep_use_pangolin ?: true,
        config.vep_flanking_seq_len ?: 10,
        config.vep_use_clinvar ?: true,
        config.vep_use_missense_zscore ?: true,
        // Stranger
        config.stranger_filter_mode ?: 'pathogenic',
        // AutoMap
        config.automap_min_roh_length ?: 500000,
        config.automap_max_roh_length ?: 100000000,
        config.automap_min_snps_in_roh ?: 20,
        config.automap_max_gap_length ?: 500000,
        true,  // automap_upd_infer_mode
        // 输出目录
        "${outdir}/qc/fastp",
        "${outdir}/qc/alignment",
        "${outdir}/qc/coverage",
        "${outdir}/alignment",
        "${outdir}/qc/sex_check",
        "${outdir}/alignment",
        "${outdir}/alignment",
        "${outdir}/alignment",
        "${outdir}/variants/snv_indel/raw",
        "${outdir}/variants/snv_indel/raw",
        "${outdir}/baf",
        "${outdir}/mt/raw",
        "${outdir}/str",
        "${outdir}/variants/snv_indel/annotated",
        "${outdir}/mt/annotated",
        "${outdir}/mei/annotated",
        "${outdir}/mei/raw",
        "${outdir}/homozygosity"
    )
}

// ============================================================================
// WES_TRIO 执行器
// ============================================================================
def run_wes_trio(config) {

    // 参考基因组通道
    def fasta = config.reference.fasta
    def fasta_fai = "${fasta}.fai"
    def fasta_dict = fasta.replaceFirst(/\.fa(sta)?$/, '.dict')

    ch_fasta = Channel.fromPath(fasta, checkIfExists: true)
    ch_fasta_fai = Channel.fromPath(fasta_fai, checkIfExists: true)
    ch_fasta_dict = Channel.fromPath(fasta_dict, checkIfExists: true)

    // BED 文件通道
    ch_target_bed = Channel.fromPath(config.target_bed, checkIfExists: true)
    ch_mt_bed = Channel.fromPath(config.mt_bed, checkIfExists: true)
    ch_snp_positions = Channel.fromPath(config.snp_positions, checkIfExists: true)

    // CNVkit 文件通道
    ch_cnvkit_reference = Channel.fromPath(config.cnvkit_reference, checkIfExists: true)
    ch_cnvkit_access_bed = config.cnvkit_access_bed ?
        Channel.fromPath(config.cnvkit_access_bed, checkIfExists: true) : Channel.empty()

    // ExpansionHunter catalog 通道
    ch_expansionhunter_catalog = Channel.fromPath(config.expansionhunter_catalog, checkIfExists: true)

    // 家系 reads 通道: [family_id, role, sample_id, [read1, read2]]
    ch_family_reads = Channel.fromList(
        config.samples.collect { sample ->
            [
                config.family_id,
                sample.role,
                sample.sample_id,
                [file(sample.read1), file(sample.read2)]
            ]
        }
    )

    // PED 文件 (如果有)
    ch_ped_file = config.ped_file ?
        Channel.fromPath(config.ped_file, checkIfExists: true) : Channel.empty()

    // 输出目录
    def outdir = config.outdir ?: params.outdir

    // 配置参数 Map
    def family_config = [
        use_gpu: config.use_gpu ?: false,
        genome_assembly: config.genome_assembly ?: 'GRCh38',
        sex_check_threshold: config.sex_check_threshold ?: 10,
        baf_max_depth: config.baf_max_depth ?: 200,
        deepvariant_model: config.deepvariant_model ?: 'WES',
        deepvariant_shards: config.deepvariant_shards,
        glnexus_config: config.glnexus_config ?: 'deepvariant',
        cnvkit_annotate: config.cnvkit_annotate ?: true,
        cnvkit_split_size: config.cnvkit_split_size ?: 5000,
        cnvkit_min_target_size: config.cnvkit_min_target_size ?: 10000,
        cnvkit_method: config.cnvkit_method ?: 'cbs',
        cnvkit_threshold: config.cnvkit_threshold,
        cnvkit_ploidy: config.cnvkit_ploidy ?: 2,
        cnvkit_drop_outliers: config.cnvkit_drop_outliers ?: true,
        expansionhunter_min_anchor: config.expansionhunter_min_anchor ?: 8,
        expansionhunter_max_irr: config.expansionhunter_max_irr ?: 100,
        tiea_min_support: config.tiea_min_support ?: 10,
        tiea_min_softclip: config.tiea_min_softclip ?: 36,
        tiea_min_mapq: config.tiea_min_mapq ?: 20,
        tiea_cluster_window: config.tiea_cluster_window ?: 10,
        tiea_threads: config.tiea_threads ?: 4,
        whatshap_chromosomes: config.whatshap_chromosomes,
        whatshap_ignore_rg: config.whatshap_ignore_rg ?: true,
        whatshap_ref_conf: config.whatshap_ref_conf ?: 20,
        vep_use_pick: config.vep_use_pick ?: true,
        vep_use_refseq_only: config.vep_use_refseq_only ?: false,
        vep_cache_dir: config.vep_cache_dir,
        vep_extra_args: config.vep_extra_args,
        vep_mei_upstream_distance: config.vep_mei_upstream_distance ?: 5000,
        vep_mei_downstream_distance: config.vep_mei_downstream_distance ?: 5000,
        vep_db_dir: config.vep_db_dir,
        vep_use_gnomad: config.vep_use_gnomad ?: true,
        vep_use_alphamissense: config.vep_use_alphamissense ?: true,
        vep_use_evoscore: config.vep_use_evoscore ?: true,
        vep_use_pangolin: config.vep_use_pangolin ?: true,
        vep_flanking_seq_len: config.vep_flanking_seq_len ?: 10,
        vep_use_clinvar: config.vep_use_clinvar ?: true,
        vep_use_missense_zscore: config.vep_use_missense_zscore ?: true,
        stranger_filter_mode: config.stranger_filter_mode ?: 'pathogenic',
        automap_min_roh_length: config.automap_min_roh_length ?: 500000,
        automap_max_roh_length: config.automap_max_roh_length ?: 100000000,
        automap_min_snps_in_roh: config.automap_min_snps_in_roh ?: 20,
        automap_max_gap_length: config.automap_max_gap_length ?: 500000,
        outdir: outdir
    ]

    // 执行 workflow
    WES_TRIO(
        ch_family_reads,
        ch_ped_file,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        ch_target_bed,
        ch_mt_bed,
        ch_snp_positions,
        ch_cnvkit_reference,
        ch_cnvkit_access_bed,
        ch_expansionhunter_catalog,
        Channel.empty(),  // vep_db_files
        family_config
    )
}

// ============================================================================
// CNV_BASELINE 执行器
// ============================================================================
def run_cnv_baseline(config) {

    // 参考基因组通道
    def fasta = config.reference.fasta
    def fasta_fai = "${fasta}.fai"

    ch_fasta = Channel.fromPath(fasta, checkIfExists: true)
    ch_fasta_fai = Channel.fromPath(fasta_fai, checkIfExists: true)

    // BWA 索引（假设在同一目录）
    def bwa_indices = Channel.fromPath(
        ["${fasta}.amb", "${fasta}.ann", "${fasta}.bwt", "${fasta}.pac", "${fasta}.sa"],
        checkIfExists: true
    )

    // BWA-MEM2 索引（可选）
    def bwamem2_indices = Channel.empty()

    // BED 文件通道
    ch_target_bed = Channel.fromPath(config.target_bed, checkIfExists: true)
    ch_access_bed = config.access_bed ?
        Channel.fromPath(config.access_bed, checkIfExists: true) : Channel.empty()

    // 样本 reads 通道: [sample_id, [read1, read2]]
    ch_samples = Channel.fromList(
        config.samples.collect { sample ->
            [
                sample.sample_id,
                [file(sample.read1), file(sample.read2)]
            ]
        }
    )

    // 判断是否使用 BWA-MEM2
    def use_bwamem2 = config.use_bwamem2 ?: false
    if (!use_bwamem2) {
        // 根据内存自动选择
        def max_mem_bytes = (params.max_memory as nextflow.util.MemoryUnit).toBytes()
        def threshold = 64L * 1024L * 1024L * 1024L  // 64GB
        use_bwamem2 = max_mem_bytes >= threshold
    }

    // CNVkit 参数
    def cnvkit = config.cnvkit ?: {}

    // 执行 workflow
    CNV_BASELINE(
        ch_samples,
        ch_fasta,
        ch_fasta_fai,
        bwa_indices,
        bwamem2_indices,
        ch_target_bed,
        config.access_bed ? true : false,
        ch_access_bed,
        use_bwamem2,
        cnvkit.annotate ?: true,
        cnvkit.split_size ?: 5000,
        cnvkit.min_target_size ?: 10000,
        cnvkit.is_female_reference ?: true
    )
}