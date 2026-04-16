#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

include { WES_SINGLE } from './workflows/wes_single'
include { QC_ALIGNMENT } from './workflows/qc_alignment'
include { CNV_BASELINE } from './workflows/cnv_baseline'

def parseConfig(configFile) {
    def jsonSlurper = new JsonSlurper()
    def config = jsonSlurper.parse(file(configFile))
    return config
}

def validateConfig(config) {
    if (!config.sample_id) {
        error "ERROR: 配置文件中未找到 sample_id"
    }
    if (!config.read1 || !config.read2) {
        error "ERROR: 配置文件中未找到 read1 或 read2"
    }
    if (!config.reference?.fasta) {
        error "ERROR: 配置文件中未找到 reference.fasta"
    }
}

// 输出目录默认路径生成函数
def getDefaultOutputDirs(outdir) {
    return [
        fastp: "${outdir}/01.QC",
        collectqc: "${outdir}/01.QC",
        xamdst: "${outdir}/01.QC",
        samtools_index: "${outdir}/02.Alignment",
        sex_check: "${outdir}/01.QC",
        markdup: "${outdir}/02.Alignment",
        bwamem: "${outdir}/02.Alignment",
        pb_fq2bam: "${outdir}/02.Alignment",
        deepvariant: "${outdir}/03.Variant/DeepVariant",
        pb_deepvariant: "${outdir}/03.Variant/DeepVariant",
        baf: "${outdir}/03.Variant/BAF",
        mutect2_mt: "${outdir}/03.Variant/Mitochondria",
        whatshap: "${outdir}/03.Variant/Phasing",
        cnvkit: "${outdir}/04.CNV",
        expansionhunter: "${outdir}/05.Annotations/STR",
        stranger: "${outdir}/05.Annotations/STR",
        vep: "${outdir}/05.Annotations/SNP_InDel",
        vep_mt: "${outdir}/05.Annotations/Mitochondria",
        vep_mei: "${outdir}/05.Annotations/MEI",
        tiea: "${outdir}/05.Other/MEI"
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
        config.output_format ?: 'cram',     // output_format
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
        outputDirs.tiea                           // tiea_output_dir
    )
}

// =============================================================================
// CNV_BASELINE - CNV 基线构建流程
// =============================================================================
workflow CNV_BASELINE {
    if (!params.config) {
        error "ERROR: 请提供配置文件 --config sample.json"
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

    // 准备参考基因组通道
    def fasta = config.reference.fasta
    ch_fasta     = Channel.value(file(fasta))
    ch_fasta_fai = Channel.value(file("${fasta}.fai"))

    // 可选文件通道
    ch_annotate = config.annotate ? Channel.value(file(config.annotate)) : Channel.empty()
    ch_access   = config.access_bed ? Channel.value(file(config.access_bed)) : Channel.empty()

    // 准备样本 BAM/CRAM 通道
    ch_alignments = Channel.create()
    config.samples.each { sample ->
        def sample_id = sample.sample_id ?: sample.id
        def alignment = sample.bam ?: sample.cram
        def index = sample.bai ?: sample.crai

        if (alignment) {
            ch_alignments = ch_alignments.mix(
                Channel.of([
                    [id: sample_id],
                    file(alignment),
                    file(index)
                ])
            )
        }
    }

    // 检查至少有一个样本
    if (ch_alignments.toList().size() == 0) {
        error "ERROR: 至少需要一个正常样本用于构建 CNV 基线"
    }

    // 获取 target BED
    ch_target_bed = Channel.value(file(config.target_bed))

    // 运行 CNV 基线构建
    CNV_BASELINE (
        ch_alignments,
        ch_fasta,
        ch_fasta_fai,
        ch_target_bed,
        ch_annotate,
        ch_access
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
