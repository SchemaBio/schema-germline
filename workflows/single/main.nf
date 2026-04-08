/*
 * WES Single Sample Pipeline
 *
 * 功能：单人全外显子测序 (WES) 分析流程
 * 输入：样本 FASTQ 文件 + 参考基因组 + 捕获区域 BED
 * 输出：变异检测结果 (SNV/Indel/CNV/STR/MEI) + QC 报告
 *
 */

// ============================================================================
// Include Modules
// ============================================================================

include { FASTP } from '../../modules/fastp/main'
include { PB_FQ2BAM } from '../../modules/parabricks/main'
include { BWAMEM; BWAMEM2 } from '../../modules/bwamem/main'
include { MARKDUPLICATES } from '../../modules/gatk/main'
include { SAMTOOLS_INDEX } from '../../modules/samtools/main'

// ============================================================================
// Workflow Definition
// ============================================================================

workflow WES_SINGLE {

    take:
    ch_reads           // Channel of [sample_id, reads] tuples
    ch_fasta           // 参考基因组 FASTA
    ch_fasta_fai       // 参考基因组索引
    val output_format  // 输出格式: 'cram' 或 'bam' (string)
    val rgid           // Read Group ID (string)
    val use_gpu        // 是否使用 GPU (boolean)

    main:
    // 预定义输出通道
    ch_alignment = Channel.empty()
    ch_alignment_index = Channel.empty()
    ch_markdup_metrics = Channel.empty()

    // =========================================================================
    // Step 1: FASTQ 质控过滤
    // =========================================================================
    FASTP(ch_reads)

    // =========================================================================
    // Step 2: 序列比对 + Step 3: 比对后处理
    // =========================================================================
    // 判断比对策略：
    //   - 有 GPU → PB_FQ2BAM (Parabricks, 无需 MarkDuplicates)
    //   - 无 GPU + 资源充足 (≥64GB) → BWAMEM2 + MARKDUPLICATES
    //   - 无 GPU + 资源不足 → BWAMEM + MARKDUPLICATES

    // 判断 process_high 资源是否满足 (内存阈值 64GB)
    mem_threshold = 64L * 1024L * 1024L * 1024L  // 64GB
    max_mem_bytes = (params.max_memory as nextflow.util.MemoryUnit).toBytes()
    has_high_resources = max_mem_bytes >= mem_threshold

    if (use_gpu) {
        // GPU 模式：使用 Parabricks (内置排序和索引，无需 MarkDuplicates)
        log.info "比对模式: Parabricks (GPU 加速)"
        PB_FQ2BAM(
            FASTP.out.clean_reads,
            ch_fasta,
            ch_fasta_fai,
            output_format,
            rgid
        )
        if (output_format == 'bam') {
            ch_alignment = ch_alignment.mix(PB_FQ2BAM.out.bam)
            ch_alignment_index = ch_alignment_index.mix(PB_FQ2BAM.out.bai)
        } else {
            ch_alignment = ch_alignment.mix(PB_FQ2BAM.out.cram)
            ch_alignment_index = ch_alignment_index.mix(PB_FQ2BAM.out.crai)
        }
    } else if (has_high_resources) {
        // 高资源模式：BWA-MEM2 + MARKDUPLICATES + SAMTOOLS_INDEX
        log.info "比对模式: BWA-MEM2 + MarkDuplicates (max_memory: ${params.max_memory})"
        BWAMEM2(
            FASTP.out.clean_reads,
            ch_fasta,
            ch_fasta_fai,
            output_format,
            rgid
        )
        // 根据输出格式选择对应的文件通道
        if (output_format == 'bam') {
            MARKDUPLICATES(
                BWAMEM2.out.bam,
                BWAMEM2.out.bai,
                ch_fasta
            )
        } else {
            MARKDUPLICATES(
                BWAMEM2.out.cram,
                BWAMEM2.out.crai,
                ch_fasta
            )
        }
        // 为 MarkDuplicates 输出创建索引
        SAMTOOLS_INDEX(MARKDUPLICATES.out.alignment)
        ch_alignment = ch_alignment.mix(SAMTOOLS_INDEX.out.alignment)
        ch_alignment_index = ch_alignment_index.mix(SAMTOOLS_INDEX.out.index)
        ch_markdup_metrics = ch_markdup_metrics.mix(MARKDUPLICATES.out.metrics)
    } else {
        // 低资源模式：BWA-MEM + MARKDUPLICATES + SAMTOOLS_INDEX
        log.info "比对模式: BWA-MEM + MarkDuplicates (max_memory: ${params.max_memory})"
        BWAMEM(
            FASTP.out.clean_reads,
            ch_fasta,
            ch_fasta_fai,
            output_format,
            rgid
        )
        if (output_format == 'bam') {
            MARKDUPLICATES(
                BWAMEM.out.bam,
                BWAMEM.out.bai,
                ch_fasta
            )
        } else {
            MARKDUPLICATES(
                BWAMEM.out.cram,
                BWAMEM.out.crai,
                ch_fasta
            )
        }
        SAMTOOLS_INDEX(MARKDUPLICATES.out.alignment)
        ch_alignment = ch_alignment.mix(SAMTOOLS_INDEX.out.alignment)
        ch_alignment_index = ch_alignment_index.mix(SAMTOOLS_INDEX.out.index)
        ch_markdup_metrics = ch_markdup_metrics.mix(MARKDUPLICATES.out.metrics)
    }

    // =========================================================================
    // Step 4: 覆盖度统计
    // =========================================================================
    // TODO

    // =========================================================================
    // Step 5: 性别检测
    // =========================================================================
    // TODO

    // =========================================================================
    // Step 6: SNV/Indel 变异检测
    // =========================================================================
    // TODO

    // =========================================================================
    // Step 7: CNV 检测
    // =========================================================================
    // TODO

    // =========================================================================
    // Step 8: STR 扩展检测
    // =========================================================================
    // TODO

    // =========================================================================
    // Step 9: MEI 检测
    // =========================================================================
    // TODO

    // =========================================================================
    // Step 10: 变异注释
    // =========================================================================
    // TODO

    // =========================================================================
    // Step 11: 变异过滤与评分
    // =========================================================================
    // TODO

    // =========================================================================
    // Emit Results
    // =========================================================================
    emit:
    clean_reads       = FASTP.out.clean_reads   // 过滤后的 FASTQ
    json_report       = FASTP.out.json_report   // JSON QC 报告
    html_report       = FASTP.out.html_report   // HTML QC 报告
    alignment         = ch_alignment            // BAM/CRAM 比对文件
    alignment_index   = ch_alignment_index      // 比对文件索引
    markdup_metrics   = ch_markdup_metrics      // MarkDuplicates 质控报告 (非GPU模式)
}