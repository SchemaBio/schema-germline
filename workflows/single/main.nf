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
include { PB_FQ2BAM; PB_DEEPVARIANT } from '../../modules/parabricks/main'
include { BWAMEM; BWAMEM2 } from '../../modules/bwamem/main'
include { MARKDUPLICATES; COLLECTQCMETRICS; MUTECT2_MT } from '../../modules/gatk/main'
include { SAMTOOLS_INDEX; SEX_CHECK_SRY } from '../../modules/samtools/main'
include { XAMDST } from '../../modules/xamdst/main'
include { BCFTOOLS_BAF_MATRIX } from '../../modules/bcftools/main'
include { DEEPVARIANT } from '../../modules/deepvariant/main'

// ============================================================================
// Workflow Definition
// ============================================================================

workflow WES_SINGLE {

    take:
    ch_reads                 // Channel of [sample_id, reads] tuples
    ch_fasta                 // 参考基因组 FASTA
    ch_fasta_fai             // 参考基因组索引 (.fai)
    ch_fasta_dict            // 参考基因组字典 (.dict)
    ch_target_bed            // 捕获区域 BED 文件 (常染色体/全外显子)
    ch_mt_bed                // 线粒体区域 BED 文件
    ch_snp_positions         // SNP 位点文件 (VCF/BED)，用于 BAF 计算
    val output_format        // 输出格式: 'cram' 或 'bam' (string)
    val rgid                 // Read Group ID (string)
    val use_gpu              // 是否使用 GPU (boolean)
    val genome_assembly      // 基因组版本: 'GRCh37' 或 'GRCh38'
    val sex_check_threshold  // 性别检测 SRY reads 阈值 (integer)
    val baf_max_depth        // BAF 计算最大深度 (integer, 默认 200)
    val deepvariant_model    // DeepVariant 模型类型: WES/WGS/PACBIO/ONT_R104 (默认 WES)
    val deepvariant_shards   // DeepVariant 并行分片数 (integer, 默认使用 task.cpus)

    main:
    // 预定义输出通道
    ch_alignment = Channel.empty()
    ch_alignment_index = Channel.empty()
    ch_markdup_metrics = Channel.empty()
    ch_vcf = Channel.empty()
    ch_vcf_tbi = Channel.empty()
    ch_gvcf = Channel.empty()
    ch_gvcf_tbi = Channel.empty()
    ch_variant_report = Channel.empty()

    // =========================================================================
    // Step 1: FASTQ 质控过滤
    // =========================================================================
    FASTP(ch_reads)

    // =========================================================================
    // Step 2: 序列比对 + Step 3: 比对后处理
    // =========================================================================
    // 判断比对策略：
    //   - 有 GPU → PB_FQ2BAM (Parabricks, 自带 MarkDup)
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
    // Step 4: 覆盖度统计 (常染色体/全外显子)
    // =========================================================================
    // 从 FASTP 输出提取 sample_id，与 alignment 组合成 XAMDST 需要的 tuple 格式
    ch_sample_id = FASTP.out.clean_reads.map { sample_id, reads -> sample_id }
    ch_alignment_tuple = ch_alignment.combine(ch_sample_id).map { alignment, sample_id -> tuple(sample_id, alignment) }

    XAMDST(ch_alignment_tuple, ch_target_bed)

    // =========================================================================
    // Step 4b: 覆盖度统计 (线粒体)
    // =========================================================================
    // 使用带 `_mt` 后缀的 sample_id 以区分线粒体统计结果
    ch_sample_id_mt = ch_sample_id.map { sample_id -> sample_id + '_mt' }
    ch_alignment_tuple_mt = ch_alignment.combine(ch_sample_id_mt).map { alignment, sample_id_mt -> tuple(sample_id_mt, alignment) }

    XAMDST.mt(ch_alignment_tuple_mt, ch_mt_bed)

    // =========================================================================
    // Step 5: 比对质量统计 (GATK CollectMultipleMetrics)
    // =========================================================================
    COLLECTQCMETRICS(
        ch_alignment,
        ch_alignment_index,
        ch_fasta,
        ch_target_bed
    )

    // =========================================================================
    // Step 6: 性别检测 (基于 SRY 区域 reads 数)
    // =========================================================================
    SEX_CHECK_SRY(
        ch_alignment,
        ch_alignment_index,
        genome_assembly,
        sex_check_threshold
    )

    // =========================================================================
    // Step 7: BAF 矩阵计算 (基于 SNP 位点)
    // =========================================================================
    BCFTOOLS_BAF_MATRIX(
        ch_alignment,
        ch_alignment_index,
        ch_fasta,
        ch_fasta_fai,
        ch_snp_positions,
        baf_max_depth
    )

    // =========================================================================
    // Step 8: SNV/Indel 变异检测 (DeepVariant)
    // =========================================================================
    // 判断变异检测策略：
    //   - 有 GPU → PB_DEEPVARIANT (Parabricks GPU 加速)
    //   - 无 GPU → DEEPVARIANT (CPU 模式)

    if (use_gpu) {
        // GPU 模式：使用 Parabricks DeepVariant
        log.info "变异检测: Parabricks DeepVariant (GPU 加速)"
        PB_DEEPVARIANT(
            ch_alignment,
            ch_alignment_index,
            ch_fasta,
            ch_fasta_fai,
            ch_fasta_dict,
            ch_target_bed,
            deepvariant_model,
            deepvariant_shards
        )
        ch_vcf = ch_vcf.mix(PB_DEEPVARIANT.out.vcf)
        ch_vcf_tbi = ch_vcf_tbi.mix(PB_DEEPVARIANT.out.vcf_tbi)
        ch_gvcf = ch_gvcf.mix(PB_DEEPVARIANT.out.gvcf)
        ch_gvcf_tbi = ch_gvcf_tbi.mix(PB_DEEPVARIANT.out.gvcf_tbi)
    } else {
        // CPU 模式：使用标准 DeepVariant
        log.info "变异检测: DeepVariant (CPU 模式)"
        DEEPVARIANT(
            ch_alignment,
            ch_alignment_index,
            ch_fasta,
            ch_fasta_fai,
            ch_fasta_dict,
            ch_target_bed,
            deepvariant_model,
            deepvariant_shards
        )
        ch_vcf = ch_vcf.mix(DEEPVARIANT.out.vcf)
        ch_vcf_tbi = ch_vcf_tbi.mix(DEEPVARIANT.out.vcf_tbi)
        ch_gvcf = ch_gvcf.mix(DEEPVARIANT.out.gvcf)
        ch_gvcf_tbi = ch_gvcf_tbi.mix(DEEPVARIANT.out.gvcf_tbi)
    }

    // =========================================================================
    // Step 8b: 线粒体变异检测 (GATK Mutect2 线粒体模式)
    // =========================================================================
    MUTECT2_MT(
        ch_alignment,
        ch_alignment_index,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        genome_assembly
    )

    // =========================================================================
    // Step 9: CNV 检测
    // =========================================================================
    // TODO

    // =========================================================================
    // Step 10: STR 扩展检测
    // =========================================================================
    // TODO

    // =========================================================================
    // Step 11: MEI 检测
    // =========================================================================
    // TODO

    // =========================================================================
    // Step 12: 变异注释
    // =========================================================================
    // TODO

    // =========================================================================
    // Step 13: 变异过滤与评分
    // =========================================================================
    // TODO

    // =========================================================================
    // Emit Results
    // =========================================================================
    emit:
    clean_reads             = FASTP.out.clean_reads              // 过滤后的 FASTQ
    json_report             = FASTP.out.json_report              // JSON QC 报告
    html_report             = FASTP.out.html_report              // HTML QC 报告
    alignment               = ch_alignment                       // BAM/CRAM 比对文件
    alignment_index         = ch_alignment_index                 // 比对文件索引
    markdup_metrics         = ch_markdup_metrics                 // MarkDuplicates 质控报告 (非GPU模式)
    coverage_report         = XAMDST.out.coverage_report         // 覆盖度统计报告 (常染色体/全外显子)
    coverage_report_mt      = XAMDST.mt.out.coverage_report      // 覆盖度统计报告 (线粒体)
    qc_metrics              = COLLECTQCMETRICS.out.metrics       // GATK QC 指标文件
    qc_pdf                  = COLLECTQCMETRICS.out.pdf           // GATK QC PDF 报告 (可选)
    sex_check_json          = SEX_CHECK_SRY.out.json             // 性别检测结果 JSON
    baf_matrix_tsv          = BCFTOOLS_BAF_MATRIX.out.baf_matrix // BAF 矩阵 TSV
    baf_matrix_json         = BCFTOOLS_BAF_MATRIX.out.baf_json   // BAF 矩阵 JSON
    vcf                     = ch_vcf                             // SNV/Indel VCF 文件 (常染色体)
    vcf_tbi                 = ch_vcf_tbi                         // VCF 索引
    gvcf                    = ch_gvcf                            // gVCF 文件
    gvcf_tbi                = ch_gvcf_tbi                        // gVCF 索引
    variant_report          = ch_variant_report                  // DeepVariant 可视化报告 (仅CPU模式)
    mt_vcf                  = MUTECT2_MT.out.vcf                 // 线粒体变异 VCF 文件
    mt_vcf_tbi              = MUTECT2_MT.out.vcf_tbi             // 线粒体 VCF 索引
    mt_stats                = MUTECT2_MT.out.stats               // 线粒体变异统计
}