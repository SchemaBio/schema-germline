#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

include { WES_SINGLE } from './workflows/wes_single'

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
