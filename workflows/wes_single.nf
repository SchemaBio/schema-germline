/*
 * WES Single Sample Pipeline
 * 
 * 流程：FASTP -> BWA_MEM2 -> GATK_MARKDUPLICATES -> CRAM
 */

// 引入模块
include { FASTP              } from '../modules/local/fastp/main'
include { BWA_MEM2           } from '../modules/local/bwa_mem2/main'
include { GATK_MARKDUPLICATES } from '../modules/local/gatk/main'

/*
 * 主工作流
 */
workflow WES_SINGLE {

    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]
    ch_bwamem2_index  // channel: path(bwamem2_index)
    ch_bwa_index      // channel: path(bwa_index)
    ch_fasta          // channel: path(fasta)
    ch_fasta_fai      // channel: path(fasta_fai)

    main:
    ch_versions = Channel.empty()

    //
    // STEP 1: FASTP - 质控过滤
    //
    FASTP ( ch_reads )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    //
    // STEP 2: BWA_MEM2 - 比对到参考基因组
    //
    BWA_MEM2 (
        FASTP.out.reads,
        ch_bwamem2_index,
        ch_bwa_index,
        ch_fasta,
        'cram'  // 输出格式
    )
    ch_versions = ch_versions.mix(BWA_MEM2.out.versions.first())

    //
    // STEP 3: GATK_MARKDUPLICATES - 标记重复
    //
    GATK_MARKDUPLICATES (
        BWA_MEM2.out.alignment,
        ch_fasta,
        'cram'  // 输出格式
    )
    ch_versions = ch_versions.mix(GATK_MARKDUPLICATES.out.versions.first())

    emit:
    reads     = FASTP.out.reads                    // channel: [ val(meta), [ reads ] ]
    fastp_json = FASTP.out.json                    // channel: [ val(meta), json ]
    fastp_html = FASTP.out.html                    // channel: [ val(meta), html ]
    cram      = GATK_MARKDUPLICATES.out.alignment  // channel: [ val(meta), cram, crai ]
    metrics   = GATK_MARKDUPLICATES.out.metrics    // channel: [ val(meta), metrics ]
    versions  = ch_versions                        // channel: versions.yml
}
