#!/usr/bin/env nextflow
/*
 * Schema Germline Pipeline - Main Entry
 */

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

// 引入工作流
include { WES_SINGLE } from './workflows/wes_single'

/*
 * 解析 JSON 配置文件
 */
def parseConfig(configFile) {
    def jsonSlurper = new JsonSlurper()
    def config = jsonSlurper.parse(file(configFile))
    return config
}

/*
 * 参数校验
 */
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
    if (!config.reference?.bwamem2_index && !config.reference?.bwa_index) {
        error "ERROR: 配置文件中未找到 reference.bwamem2_index 或 reference.bwa_index"
    }
}

/*
 * 主流程
 */
workflow {
    // 解析配置文件
    if (!params.config) {
        error "ERROR: 请提供配置文件 --config sample.json"
    }
    
    def config = parseConfig(params.config)
    validateConfig(config)

    // 设置输出目录
    params.outdir = config.outdir ?: params.outdir

    // 构建样本 channel
    def meta = [
        id: config.sample_id,
        read_group: "@RG\\tID:${config.sample_id}\\tSM:${config.sample_id}\\tPL:ILLUMINA"
    ]
    def reads = [file(config.read1), file(config.read2)]
    ch_reads = Channel.of([meta, reads])
    
    // 参考基因组 channels
    def ref = config.reference
    ch_fasta         = Channel.value(file(ref.fasta))
    ch_fasta_fai     = ref.fasta_fai ? Channel.value(file(ref.fasta_fai)) : Channel.value(file("${ref.fasta}.fai"))
    ch_bwamem2_index = ref.bwamem2_index ? Channel.value(file(ref.bwamem2_index)) : Channel.value([])
    ch_bwa_index     = ref.bwa_index ? Channel.value(file(ref.bwa_index)) : Channel.value([])

    // 运行 WES 单样本流程
    WES_SINGLE (
        ch_reads,
        ch_bwamem2_index,
        ch_bwa_index,
        ch_fasta,
        ch_fasta_fai
    )
}
