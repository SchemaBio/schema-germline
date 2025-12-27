#!/usr/bin/env nextflow
/*
 * Schema Germline Pipeline - Main Entry
 */

nextflow.enable.dsl = 2

// 引入工作流
include { WES_SINGLE } from './workflows/wes_single'

/*
 * 参数校验
 */
def validateParams() {
    if (!params.input) {
        error "ERROR: 请提供输入样本表 --input samplesheet.csv"
    }
    if (!params.fasta) {
        error "ERROR: 请提供参考基因组 --fasta reference.fa"
    }
    if (!params.bwamem2_index && !params.bwa_index) {
        error "ERROR: 请提供 BWA 索引 --bwamem2_index 或 --bwa_index"
    }
}

/*
 * 解析样本表
 * 格式: sample_id,read1,read2
 */
def parseSamplesheet(samplesheet) {
    return Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample_id,
                read_group: "@RG\\tID:${row.sample_id}\\tSM:${row.sample_id}\\tPL:ILLUMINA"
            ]
            def reads = [file(row.read1), file(row.read2)]
            return [meta, reads]
        }
}

/*
 * 主流程
 */
workflow {
    // 参数校验
    validateParams()

    // 准备输入 channels
    ch_reads = parseSamplesheet(params.input)
    
    ch_fasta         = Channel.value(file(params.fasta))
    ch_fasta_fai     = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : Channel.value(file("${params.fasta}.fai"))
    ch_bwamem2_index = params.bwamem2_index ? Channel.value(file(params.bwamem2_index)) : Channel.value([])
    ch_bwa_index     = params.bwa_index ? Channel.value(file(params.bwa_index)) : Channel.value([])

    // 运行 WES 单样本流程
    WES_SINGLE (
        ch_reads,
        ch_bwamem2_index,
        ch_bwa_index,
        ch_fasta,
        ch_fasta_fai
    )
}
