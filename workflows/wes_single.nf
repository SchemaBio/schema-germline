/*
 * WES Single Sample Pipeline
 */

include { FASTP              } from '../modules/local/fastp/main'
include { BWA_MEM            } from '../modules/local/bwa_mem2/main'
include { BWA_MEM2           } from '../modules/local/bwa_mem2/main'
include { GATK_MARKDUPLICATES } from '../modules/local/gatk/main'

workflow WES_SINGLE {

    take:
    ch_reads
    ch_fasta
    ch_fasta_fai
    // bwa 索引
    ch_bwa_amb
    ch_bwa_ann
    ch_bwa_bwt
    ch_bwa_pac
    ch_bwa_sa
    // bwa-mem2 索引
    ch_bwamem2_0123
    ch_bwamem2_bwt2bit64
    // 使用哪个比对器
    use_bwamem2

    main:
    ch_versions = Channel.empty()

    FASTP ( ch_reads )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    if (use_bwamem2) {
        BWA_MEM2 (
            FASTP.out.reads,
            ch_fasta,
            ch_fasta_fai,
            ch_bwamem2_0123,
            ch_bwa_amb,
            ch_bwa_ann,
            ch_bwamem2_bwt2bit64,
            ch_bwa_pac,
            'cram'
        )
        ch_alignment = BWA_MEM2.out.alignment
        ch_versions = ch_versions.mix(BWA_MEM2.out.versions.first())
    } else {
        BWA_MEM (
            FASTP.out.reads,
            ch_fasta,
            ch_fasta_fai,
            ch_bwa_amb,
            ch_bwa_ann,
            ch_bwa_bwt,
            ch_bwa_pac,
            ch_bwa_sa,
            'cram'
        )
        ch_alignment = BWA_MEM.out.alignment
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    }

    GATK_MARKDUPLICATES (
        ch_alignment,
        ch_fasta,
        'cram'
    )
    ch_versions = ch_versions.mix(GATK_MARKDUPLICATES.out.versions.first())

    emit:
    reads      = FASTP.out.reads
    fastp_json = FASTP.out.json
    fastp_html = FASTP.out.html
    cram       = GATK_MARKDUPLICATES.out.alignment
    metrics    = GATK_MARKDUPLICATES.out.metrics
    versions   = ch_versions
}
