#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

include { WES_SINGLE } from './workflows/wes_single'
include { WES_TRIO_PIPELINE } from './workflows/trio'
include { QC_ALIGNMENT } from './workflows/qc_alignment'
include { CNV_BASELINE } from './workflows/cnv_baseline'

def parseConfig(configFile) {
    def jsonSlurper = new JsonSlurper()
    def config = jsonSlurper.parse(file(configFile))
    return config
}

def validateConfig(config) {
    if (!config.sample_id && !config.samples) {
        error "ERROR: 配置文件中未找到 sample_id (单人) 或 samples (家系)"
    }
    if (!config.read1 || !config.read2) {
        error "ERROR: 配置文件中未找到 read1 或 read2"
    }
    if (!config.reference?.fasta) {
        error "ERROR: 配置文件中未找到 reference.fasta"
    }
}

// 输出目录默认路径生成函数
// 参考 Illumina DRAGEN / Sentieon / GATK Best Practices 的目录结构设计
// 按数据类型分类，中间文件和最终结果分离，QC报告集中管理
def getDefaultOutputDirs(outdir) {
    return [
        // QC 质控报告 (集中管理)
        fastp: "${outdir}/qc/fastp",
        collectqc: "${outdir}/qc/alignment",
        xamdst: "${outdir}/qc/coverage",
        sex_check: "${outdir}/qc/sex_check",
        // 比对结果
        samtools_index: "${outdir}/alignment",
        markdup: "${outdir}/alignment",
        bwamem: "${outdir}/alignment",
        pb_fq2bam: "${outdir}/alignment",
        // SNV/Indel 变异检测
        deepvariant: "${outdir}/variants/snv_indel/raw",
        pb_deepvariant: "${outdir}/variants/snv_indel/raw",
        whatshap: "${outdir}/variants/snv_indel/phased",
        vep: "${outdir}/variants/snv_indel/annotated",
        // BAF 矩阵 (辅助分析)
        baf: "${outdir}/baf",
        // 线粒体变异
        mutect2_mt: "${outdir}/variants/mt/raw",
        vep_mt: "${outdir}/variants/mt/annotated",
        // CNV 检测
        cnvkit: "${outdir}/variants/cnv",
        // STR 扩展检测
        expansionhunter: "${outdir}/variants/str/raw",
        stranger: "${outdir}/variants/str/annotated",
        // MEI 检测
        tiea: "${outdir}/variants/mei/raw",
        vep_mei: "${outdir}/variants/mei/annotated",
        // ROH/UPD 分析
        automap: "${outdir}/homozygosity"
    ]
}

workflow {
    if (!params.config) {
        error "ERROR: 请提供配置文件 --config sample.json"
    }
    
    def config = parseConfig(params.config)
    validateConfig(config)

    params.outdir = config.outdir ?: params.outdir

    def meta = [
        id: config.sample_id,
        read_group: "@RG\\tID:${config.sample_id}\\tSM:${config.sample_id}\\tPL:ILLUMINA"
    ]
    def reads = [file(config.read1), file(config.read2)]
    ch_reads = Channel.of([meta, reads])
    
    def fasta = config.reference.fasta
    ch_fasta     = Channel.value(file(fasta))
    ch_fasta_fai = Channel.value(file("${fasta}.fai"))
    
    // bwa 索引
    ch_bwa_amb = Channel.value(file("${fasta}.amb"))
    ch_bwa_ann = Channel.value(file("${fasta}.ann"))
    ch_bwa_bwt = Channel.value(file("${fasta}.bwt"))
    ch_bwa_pac = Channel.value(file("${fasta}.pac"))
    ch_bwa_sa  = Channel.value(file("${fasta}.sa"))
    
    // bwa-mem2 索引
    ch_bwamem2_0123      = Channel.value(file("${fasta}.0123"))
    ch_bwamem2_bwt2bit64 = Channel.value(file("${fasta}.bwt.2bit.64"))

    // 根据 max_memory 决定使用哪个比对器
    def mem_threshold = 64L * 1024L * 1024L * 1024L  // 64GB
    def max_mem_bytes = (params.max_memory as nextflow.util.MemoryUnit).toBytes()
    def use_bwamem2 = max_mem_bytes >= mem_threshold

    log.info "比对器选择: ${use_bwamem2 ? 'bwa-mem2' : 'bwa'} (max_memory: ${params.max_memory})"

    // 获取默认输出目录
    def outdir = config.outdir ?: params.outdir
    def outputDirs = getDefaultOutputDirs(outdir)

    // 准备其他必要参数的通道和值
    // 注意: 以下参数需要从配置文件或其他来源获取
    ch_target_bed = config.target_bed ? Channel.value(file(config.target_bed)) : Channel.empty()
    ch_mt_bed = config.mt_bed ? Channel.value(file(config.mt_bed)) : Channel.empty()
    ch_snp_positions = config.snp_positions ? Channel.value(file(config.snp_positions)) : Channel.empty()
    ch_fasta_dict = Channel.value(file("${fasta.replaceAll(/\\.fai$/, '')}.dict"))

    // CNVkit 参数
    ch_cnvkit_reference = config.cnvkit_reference ? Channel.value(file(config.cnvkit_reference)) : Channel.empty()
    ch_cnvkit_access_bed = config.cnvkit_access_bed ? Channel.value(file(config.cnvkit_access_bed)) : Channel.empty()

    // ExpansionHunter 参数
    ch_expansionhunter_catalog = config.expansionhunter_catalog ? Channel.value(file(config.expansionhunter_catalog)) : Channel.empty()

    WES_SINGLE (
        ch_reads,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        ch_target_bed,
        ch_mt_bed,
        ch_snp_positions,
        config.sample_id,                     // rgid
        config.use_gpu ?: false,              // use_gpu
        config.genome_assembly ?: 'GRCh38',   // genome_assembly
        config.sex_check_threshold ?: 10,     // sex_check_threshold
        config.baf_max_depth ?: 200,          // baf_max_depth
        config.deepvariant_model ?: 'WES',    // deepvariant_model
        config.deepvariant_shards ?: null,    // deepvariant_shards
        // CNVkit 参数
        ch_cnvkit_reference,
        ch_cnvkit_access_bed,
        config.cnvkit_annotate ?: true,       // cnvkit_annotate
        config.cnvkit_split_size ?: 5000,     // cnvkit_split_size
        config.cnvkit_min_target_size ?: 10000, // cnvkit_min_target_size
        config.cnvkit_method ?: 'cbs',        // cnvkit_method
        config.cnvkit_threshold ?: null,      // cnvkit_threshold
        config.cnvkit_ploidy ?: 2,            // cnvkit_ploidy
        config.cnvkit_drop_outliers ?: true,  // cnvkit_drop_outliers
        outputDirs.cnvkit,                    // cnvkit_output_dir
        // ExpansionHunter 参数
        ch_expansionhunter_catalog,
        config.expansionhunter_min_anchor ?: 8,   // expansionhunter_min_anchor
        config.expansionhunter_max_irr ?: 100,    // expansionhunter_max_irr
        outputDirs.expansionhunter,               // expansionhunter_output_dir
        // TIEA-WES 参数
        config.tiea_min_support ?: 10,            // tiea_min_support
        config.tiea_min_softclip ?: 36,           // tiea_min_softclip
        config.tiea_min_mapq ?: 20,               // tiea_min_mapq
        config.tiea_cluster_window ?: 10,         // tiea_cluster_window
        config.tiea_threads ?: 4,                 // tiea_threads
        // WhatsHap 参数
        config.whatshap_chromosomes ?: null,      // whatshap_chromosomes
        config.whatshap_ignore_rg ?: true,        // whatshap_ignore_rg
        config.whatshap_ref_conf ?: 20,           // whatshap_ref_conf
        outputDirs.whatshap,                      // whatshap_output_dir
        // VEP 参数
        config.vep_use_pick ?: true,              // vep_use_pick
        config.vep_use_refseq_only ?: false,      // vep_use_refseq_only
        config.vep_cache_dir ?: null,             // vep_cache_dir
        config.vep_extra_args ?: null,            // vep_extra_args
        config.vep_mei_upstream_distance ?: 5000, // vep_mei_upstream_distance
        config.vep_mei_downstream_distance ?: 5000, // vep_mei_downstream_distance
        // VEP 自定义数据库参数
        config.vep_db_dir ?: null,                // vep_db_dir
        config.vep_use_gnomad ?: true,            // vep_use_gnomad
        config.vep_use_alphamissense ?: true,     // vep_use_alphamissense
        config.vep_use_evoscore ?: true,          // vep_use_evoscore
        config.vep_use_civic ?: true,             // vep_use_civic
        config.vep_use_mskcc ?: true,             // vep_use_mskcc
        config.vep_use_pangolin ?: true,          // vep_use_pangolin
        // Stranger 参数
        config.stranger_filter_mode ?: 'pathogenic', // stranger_filter_mode
        // AutoMap ROH/UPD 参数
        params.automap_min_roh_length ?: 500000,     // automap_min_roh_length
        params.automap_max_roh_length ?: 100000000,  // automap_max_roh_length
        params.automap_min_snps_in_roh ?: 20,        // automap_min_snps_in_roh
        params.automap_max_gap_length ?: 500000,     // automap_max_gap_length
        true,                                        // automap_upd_infer_mode (单人模式默认启用)
        // 输出目录参数
        outputDirs.fastp,                         // fastp_output_dir
        outputDirs.collectqc,                     // collectqc_output_dir
        outputDirs.xamdst,                        // xamdst_output_dir
        outputDirs.samtools_index,                // samtools_index_output_dir
        outputDirs.sex_check,                     // sex_check_output_dir
        outputDirs.markdup,                       // markdup_output_dir
        outputDirs.bwamem,                        // bwamem_output_dir
        outputDirs.pb_fq2bam,                     // pb_fq2bam_output_dir
        outputDirs.deepvariant,                   // deepvariant_output_dir
        outputDirs.pb_deepvariant,                // pb_deepvariant_output_dir
        outputDirs.baf,                           // baf_output_dir
        outputDirs.mutect2_mt,                    // mutect2_mt_output_dir
        outputDirs.stranger,                      // stranger_output_dir
        outputDirs.vep,                           // vep_output_dir
        outputDirs.vep_mt,                        // vep_mt_output_dir
        outputDirs.vep_mei,                       // vep_mei_output_dir
        outputDirs.tiea,                           // tiea_output_dir
        outputDirs.automap                          // automap_output_dir
    )
}

// =============================================================================
// CNV_BASELINE - CNV 基线构建流程 (从 FASTQ)
// =============================================================================
workflow CNV_BASELINE {
    if (!params.config) {
        error "ERROR: 请提供配置文件 --config baseline.json"
    }

    def config = parseConfig(params.config)
    if (!config.samples) {
        error "ERROR: 配置文件中未找到 samples 数组 (用于 CNV 基线构建)"
    }
    if (!config.reference?.fasta) {
        error "ERROR: 配置文件中未找到 reference.fasta"
    }
    if (!config.target_bed) {
        error "ERROR: 配置文件中未找到 target_bed (捕获区域 BED 文件)"
    }

    params.outdir = config.outdir ?: params.outdir

    // =========================================================================
    // 准备参考基因组通道
    // =========================================================================
    def fasta = config.reference.fasta
    ch_fasta     = Channel.value(file(fasta))
    ch_fasta_fai = Channel.value(file("${fasta}.fai"))

    // BWA 索引文件集合
    ch_bwa_indices = Channel.fromPath("${fasta}.{amb,ann,bwt,pac,sa}")

    // BWA-MEM2 索引文件集合 (可选)
    ch_bwamem2_indices = Channel.fromPath("${fasta}.0123") +
                         Channel.fromPath("${fasta}.bwt.2bit.64")

    // =========================================================================
    // 准备样本 FASTQ 文件通道
    // =========================================================================
    ch_samples = Channel.create()
    config.samples.each { sample ->
        def sample_id = sample.sample_id ?: sample.id
        def read1 = sample.read1
        def read2 = sample.read2

        if (read1 && read2) {
            def meta = [
                id: sample_id,
                read_group: "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA"
            ]
            ch_samples = ch_samples.mix(
                Channel.of([meta, [file(read1), file(read2)]])
            )
        }
    }

    // 检查至少有一个样本
    if (ch_samples.toList().size() == 0) {
        error "ERROR: 至少需要一个正常样本用于构建 CNV 基线"
    }

    // =========================================================================
    // 准备其他参数
    // =========================================================================
    ch_target_bed = Channel.value(file(config.target_bed))
    def has_access_bed = config.access_bed ? true : false
    ch_access_bed = config.access_bed ? Channel.value(file(config.access_bed)) : Channel.empty()

    // 根据 max_memory 决定使用哪个比对器
    def mem_threshold = 64L * 1024L * 1024L * 1024L  // 64GB
    def max_mem_bytes = (params.max_memory as nextflow.util.MemoryUnit).toBytes()
    def use_bwamem2 = config.use_bwamem2 ?: (max_mem_bytes >= mem_threshold)

    // CNVkit 参数
    def cnvkit_config = config.cnvkit ?: [:]
    def annotate = cnvkit_config.annotate ?: true
    def split_size = cnvkit_config.split_size ?: 5000
    def min_target_size = cnvkit_config.min_target_size ?: 10000
    def is_female_reference = cnvkit_config.is_female_reference ?: true

    log.info "CNV 基线构建: ${config.samples.size()} 个样本"
    log.info "比对器选择: ${use_bwamem2 ? 'bwa-mem2' : 'bwa'}"

    // =========================================================================
    // 运行 CNV 基线构建流程
    // =========================================================================
    CNV_BASELINE (
        ch_samples,
        ch_fasta,
        ch_fasta_fai,
        ch_bwa_indices,
        ch_bwamem2_indices,
        ch_target_bed,
        has_access_bed,
        ch_access_bed,
        use_bwamem2,
        annotate,
        split_size,
        min_target_size,
        is_female_reference
    )
}

// =============================================================================
// QC_ALIGNMENT - 质控和比对流程 (简化版)
// =============================================================================
workflow QC_ALIGNMENT {
    if (!params.config) {
        error "ERROR: 请提供配置文件 --config sample.json"
    }

    def config = parseConfig(params.config)
    validateConfig(config)

    params.outdir = config.outdir ?: params.outdir

    def meta = [
        id: config.sample_id,
        read_group: "@RG\\tID:${config.sample_id}\\tSM:${config.sample_id}\\tPL:ILLUMINA"
    ]
    def reads = [file(config.read1), file(config.read2)]
    ch_reads = Channel.of([meta, reads])

    def fasta = config.reference.fasta
    ch_fasta     = Channel.value(file(fasta))
    ch_fasta_fai = Channel.value(file("${fasta}.fai"))

    // bwa 索引
    ch_bwa_amb = Channel.value(file("${fasta}.amb"))
    ch_bwa_ann = Channel.value(file("${fasta}.ann"))
    ch_bwa_bwt = Channel.value(file("${fasta}.bwt"))
    ch_bwa_pac = Channel.value(file("${fasta}.pac"))
    ch_bwa_sa  = Channel.value(file("${fasta}.sa"))

    // bwa-mem2 索引
    ch_bwamem2_0123      = Channel.value(file("${fasta}.0123"))
    ch_bwamem2_bwt2bit64 = Channel.value(file("${fasta}.bwt.2bit.64"))

    // 根据 max_memory 决定使用哪个比对器
    def mem_threshold = 64L * 1024L * 1024L * 1024L  // 64GB
    def max_mem_bytes = (params.max_memory as nextflow.util.MemoryUnit).toBytes()
    def use_bwamem2 = max_mem_bytes >= mem_threshold

    log.info "比对器选择: ${use_bwamem2 ? 'bwa-mem2' : 'bwa'} (max_memory: ${params.max_memory})"

    QC_ALIGNMENT (
        ch_reads,
        ch_fasta,
        ch_fasta_fai,
        ch_bwa_amb,
        ch_bwa_ann,
        ch_bwa_bwt,
        ch_bwa_pac,
        ch_bwa_sa,
        ch_bwamem2_0123,
        ch_bwamem2_bwt2bit64,
        use_bwamem2
    )
}

// =============================================================================
// WES_TRIO - 家系全外显子测序分析流程
// =============================================================================
workflow WES_TRIO {
    if (!params.config) {
        error "ERROR: 请提供配置文件 --config trio.json"
    }

    def config = parseConfig(params.config)

    // 验证家系配置
    if (!config.samples || config.samples.size() < 3) {
        error "ERROR: 家系流程需要至少3个样本 (proband + father + mother)"
    }
    if (!config.family_id) {
        error "ERROR: 配置文件中未找到 family_id"
    }
    if (!config.reference?.fasta) {
        error "ERROR: 配置文件中未找到 reference.fasta"
    }

    // 验证家系成员角色
    def roles = config.samples.collect { it.role }
    if (!roles.contains('proband')) {
        error "ERROR: 家系配置中未找到 proband (先证者)"
    }
    if (!roles.contains('father')) {
        error "ERROR: 家系配置中未找到 father (父)"
    }
    if (!roles.contains('mother')) {
        error "ERROR: 家系配置中未找到 mother (母)"
    }

    params.outdir = config.outdir ?: params.outdir

    def family_id = config.family_id
    log.info "家系分析: ${family_id}"

    // 准备参考基因组通道
    def fasta = config.reference.fasta
    ch_fasta     = Channel.value(file(fasta))
    ch_fasta_fai = Channel.value(file("${fasta}.fai"))
    ch_fasta_dict = Channel.value(file("${fasta.replaceAll(/\\.fa$/, '')}.dict"))

    // 准备捕获区域通道
    ch_target_bed = config.target_bed ? Channel.value(file(config.target_bed)) : Channel.empty()
    ch_mt_bed = config.mt_bed ? Channel.value(file(config.mt_bed)) : Channel.empty()
    ch_snp_positions = config.snp_positions ? Channel.value(file(config.snp_positions)) : Channel.empty()

    // CNVkit 参数
    ch_cnvkit_reference = config.cnvkit_reference ? Channel.value(file(config.cnvkit_reference)) : Channel.empty()
    ch_cnvkit_access_bed = config.cnvkit_access_bed ? Channel.value(file(config.cnvkit_access_bed)) : Channel.empty()

    // ExpansionHunter 参数
    ch_expansionhunter_catalog = config.expansionhunter_catalog ? Channel.value(file(config.expansionhunter_catalog)) : Channel.empty()

    // 准备家系成员 reads 通道
    // 格式: [family_id, member_role, sample_id, [read1, read2]]
    ch_family_reads = Channel.create()
    config.samples.each { sample ->
        def sample_id = sample.sample_id
        def role = sample.role
        def reads = [file(sample.read1), file(sample.read2)]

        ch_family_reads = ch_family_reads.mix(
            Channel.of([family_id, role, sample_id, reads])
        )
    }

    // 创建 PED 文件
    def ped_file = createPedFile(config)

    // 获取输出目录
    def outdir = config.outdir ?: params.outdir
    def outputDirs = getDefaultTrioOutputDirs(outdir, family_id)

    // 判断比对策略
    def mem_threshold = 64L * 1024L * 1024L * 1024L  // 64GB
    def max_mem_bytes = (params.max_memory as nextflow.util.MemoryUnit).toBytes()
    def has_high_resources = max_mem_bytes >= mem_threshold
    def use_gpu = config.use_gpu ?: false

    log.info "比对模式: ${use_gpu ? 'Parabricks (GPU)' : (has_high_resources ? 'BWA-MEM2' : 'BWA-MEM')}"

    // 构建配置参数 Map
    def trio_config = [
        use_gpu: use_gpu,
        genome_assembly: config.genome_assembly ?: 'GRCh38',
        sex_check_threshold: config.sex_check_threshold ?: 10,
        baf_max_depth: config.baf_max_depth ?: 200,
        deepvariant_model: config.deepvariant_model ?: 'WES',
        deepvariant_shards: config.deepvariant_shards,
        glnexus_config: params.glnexus_config ?: 'deepvariant',
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
        automap_min_roh_length: params.automap_min_roh_length ?: 500000,
        automap_max_roh_length: params.automap_max_roh_length ?: 100000000,
        automap_min_snps_in_roh: params.automap_min_snps_in_roh ?: 20,
        automap_max_gap_length: params.automap_max_gap_length ?: 500000,
        outdir: "${outdir}/${family_id}"
    ]

    WES_TRIO_PIPELINE (
        ch_family_reads,
        Channel.value(ped_file),
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        ch_target_bed,
        ch_mt_bed,
        ch_snp_positions,
        ch_cnvkit_reference,
        ch_cnvkit_access_bed,
        ch_expansionhunter_catalog,
        Channel.empty(),        // ch_vep_db_files (可选)
        trio_config
    )
}

// 创建 PED 文件的辅助函数
def createPedFile(config) {
    def family_id = config.family_id
    def ped_path = "${config.outdir ?: params.outdir}/${family_id}/family.ped"

    // 解析家系成员信息
    def proband = config.samples.find { it.role == 'proband' }
    def father = config.samples.find { it.role == 'father' }
    def mother = config.samples.find { it.role == 'mother' }

    // PED 格式: FamilyID SampleID FatherID MotherID Sex Phenotype
    // Sex: 1=male, 2=female, 0=unknown
    // Phenotype: 1=unaffected, 2=affected, 0=unknown
    def ped_content = """
${family_id} ${proband.sample_id} ${father.sample_id} ${mother.sample_id} 0 2
${family_id} ${father.sample_id} 0 0 1 1
${family_id} ${mother.sample_id} 0 0 2 1
"""

    // 写入 PED 文件
    def ped_file = new File(ped_path)
    ped_file.parentFile.mkdirs()
    ped_file.text = ped_content.trim()

    return ped_file
}

// Trio 输出目录生成函数
def getDefaultTrioOutputDirs(outdir, family_id) {
    def base = "${outdir}/${family_id}"
    return [
        qc: "${base}/qc",
        fastp: "${base}/qc/fastp",
        alignment_qc: "${base}/qc/alignment",
        coverage: "${base}/qc/coverage",
        sex_check: "${base}/qc/sex_check",
        peddy: "${base}/qc/peddy",
        alignment: "${base}/alignment",
        snv_indel_raw: "${base}/variants/snv_indel/raw",
        snv_indel_joint: "${base}/variants/snv_indel/joint",
        snv_indel_phased: "${base}/variants/snv_indel/phased",
        snv_indel_annotated: "${base}/variants/snv_indel/annotated",
        de_novo: "${base}/variants/de_novo",
        cnv: "${base}/proband/cnv",
        str_raw: "${base}/proband/str/raw",
        str_annotated: "${base}/proband/str/annotated",
        str_filtered: "${base}/proband/str/filtered",
        mei_raw: "${base}/proband/mei/raw",
        mei_annotated: "${base}/proband/mei/annotated",
        mt_raw: "${base}/proband/mt/raw",
        mt_annotated: "${base}/proband/mt/annotated",
        baf: "${base}/proband/baf",
        homozygosity: "${base}/homozygosity"
    ]
}
