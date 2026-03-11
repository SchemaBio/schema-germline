/*
 * QC & Alignment Pipeline - Simplified Workflow
 *
 * 功能：只进行质控和比对，输出 CRAM 文件
 *
 * 输出目录结构:
 * ${params.outdir}/
 * └── {sample_id}/
 *     ├── 01_qc/              # 质控报告 (fastp)
 *     └── 02_alignment/       # 比对结果 (cram files)
 *
 * 使用方式:
 *   nextflow run main.nf -profile docker -entry QC_ALIGNMENT --config sample.json
 */

// ============================================================================
// Include Modules
// ============================================================================

// Quality Control
include { FASTP } from '../modules/local/fastp/main'

// Alignment
include { BWA_MEM } from '../modules/local/bwa_mem2/main'
include { BWA_MEM2 } from '../modules/local/bwa_mem2/main'


// ============================================================================
// Workflow Definition
// ============================================================================

workflow QC_ALIGNMENT {

    take:
    ch_reads                       // [meta, [read1, read2]]
    ch_fasta                       // 参考基因组
    ch_fasta_fai                   // 参考基因组索引
    // BWA 索引
    ch_bwa_amb
    ch_bwa_ann
    ch_bwa_bwt
    ch_bwa_pac
    ch_bwa_sa
    // BWA-MEM2 索引
    ch_bwamem2_0123
    ch_bwamem2_bwt2bit64
    // 是否使用 BWA-MEM2
    use_bwamem2

    main:
    ch_versions = Channel.empty()

    // =========================================================================
    // Step 1: FASTP - Quality Control & Trimming
    // =========================================================================
    FASTP(ch_reads)
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    // =========================================================================
    // Step 2: BWA - Alignment
    // =========================================================================
    if (use_bwamem2) {
        BWA_MEM2(
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
        BWA_MEM(
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

    // =========================================================================
    // Emit Final Results
    // =========================================================================
    emit:
    // Sample info
    sample_id = FASTP.out.reads.map { meta, reads -> meta.id }

    // Raw data processing
    reads = FASTP.out.reads
    fastp_json = FASTP.out.json
    fastp_html = FASTP.out.html

    // Alignment (sorted.cram)
    cram = ch_alignment

    // Versions
    versions = ch_versions
}