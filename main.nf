#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

include { WES_SINGLE } from './workflows/wes_single'
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

    WES_SINGLE (
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
