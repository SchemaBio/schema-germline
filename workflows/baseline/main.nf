/*
 * CNV Baseline Pipeline - 从 FASTQ 构建 CNVkit 参考基线
 *
 * 功能：
 *   1. 输入多个正常样本的 FASTQ 文件
 *   2. 使用 BWA-MEM 或 BWA-MEM2 进行比对
 *   3. 所有样本比对完成后，使用 CNVkit 计算覆盖度
 *   4. 合并所有覆盖度文件生成 reference.cnn 基线
 *
 * 输入：JSON 配置文件，包含多个样本信息
 * 输出：reference.cnn 基线文件
 *
 * 使用方法：
 *   nextflow run main.nf -entry CNV_BASELINE \
 *       --config baseline.json \
 *       -profile docker
 *
 * 配置文件格式 (JSON)：
 *   {
 *     "outdir": "./results/baseline",
 *     "reference": {
 *       "fasta": "/path/to/reference.fa"
 *     },
 *     "target_bed": "/path/to/capture.bed",
 *     "access_bed": "/path/to/access.bed",     // 可选
 *     "use_bwamem2": true,                      // 可选，默认根据内存自动选择
 *     "cnvkit": {
 *       "annotate": true,
 *       "split_size": 5000,
 *       "min_target_size": 10000,
 *       "is_female_reference": true
 *     },
 *     "samples": [
 *       {
 *         "sample_id": "normal1",
 *         "read1": "/path/to/normal1_R1.fastq.gz",
 *         "read2": "/path/to/normal1_R2.fastq.gz"
 *       },
 *       {
 *         "sample_id": "normal2",
 *         "read1": "/path/to/normal2_R1.fastq.gz",
 *         "read2": "/path/to/normal2_R2.fastq.gz"
 *       }
 *     ]
 *   }
 *
 * 说明：
 *   - 至少需要 1 个正常样本，推荐 >= 5 个以获得稳定基线
 *   - 支持 BWA-MEM 和 BWA-MEM2 (根据内存自动选择)
 *   - 可选使用 Fastp 进行质控和接头去除
 */

nextflow.enable.dsl = 2

// ============================================================================
// Include Modules
// ============================================================================

include { BWAMEM } from '../modules/bwamem'
include { BWAMEM2 } from '../modules/bwamem'
include { SAMTOOLS_INDEX } from '../modules/samtools'
include { CNVKIT_TARGET_ANTITARGET } from '../modules/cnvkit'
include { CNVKIT_COVERAGE } from '../modules/cnvkit'
include { CNVKIT_REFERENCE_BUILD } from '../modules/cnvkit'

// ============================================================================
// Workflow Definition
// ============================================================================

workflow CNV_BASELINE {

    take:
    ch_samples           // Channel of [meta, [read1, read2]] tuples
    ch_fasta             // 参考基因组 FASTA
    ch_fasta_fai         // 参考基因组索引 (.fai)
    ch_bwa_indices       // BWA 索引文件集合
    ch_bwamem2_indices   // BWA-MEM2 索引文件集合 (可选)
    ch_target_bed        // 用户捕获区域 BED
    val has_access_bed   // 是否有 access_bed 文件
    ch_access_bed        // 可访问区域 BED (可选)
    val use_bwamem2      // 是否使用 BWA-MEM2
    val annotate         // 是否添加基因注释
    val split_size       // 分割大小
    val min_target_size  // 最小目标大小
    val is_female_reference // 是否使用女性参考

    main:
    // =========================================================================
    // Step 1: 比对 - BWA-MEM 或 BWA-MEM2
    // =========================================================================
    if (use_bwamem2) {
        BWAMEM2(
            ch_samples,
            ch_fasta,
            ch_bwamem2_indices
        )
        ch_bam_raw = BWAMEM2.out.bam.map { meta, bam -> bam }
    } else {
        BWAMEM(
            ch_samples,
            ch_fasta,
            ch_bwa_indices
        )
        ch_bam_raw = BWAMEM.out.bam.map { meta, bam -> bam }
    }

    // =========================================================================
    // Step 2: 创建 BAM 索引
    // =========================================================================
    SAMTOOLS_INDEX(ch_bam_raw)

    // =========================================================================
    // Step 3: 生成 target/antitarget BED
    // =========================================================================
    // 处理可选的 access_bed：如果没有 access_bed，使用占位符
    ch_access_bed_input = has_access_bed ? ch_access_bed : Channel.value('NO_ACCESS_BED')

    CNVKIT_TARGET_ANTITARGET(
        ch_target_bed,
        has_access_bed,
        ch_access_bed_input,
        annotate,
        split_size,
        min_target_size
    )

    // =========================================================================
    // Step 4: 计算每个样本的覆盖度 (并行)
    // =========================================================================
    CNVKIT_COVERAGE(
        SAMTOOLS_INDEX.out.alignment,
        SAMTOOLS_INDEX.out.index,
        CNVKIT_TARGET_ANTITARGET.out.target_bed,
        CNVKIT_TARGET_ANTITARGET.out.antitarget_bed
    )

    // =========================================================================
    // Step 5: 收集所有覆盖度文件并构建基线
    // =========================================================================
    // collect() 将所有元素收集为一个列表，等待所有样本完成
    ch_target_coverages = CNVKIT_COVERAGE.out.target_coverage.collect()
    ch_antitarget_coverages = CNVKIT_COVERAGE.out.antitarget_coverage.collect()

    CNVKIT_REFERENCE_BUILD(
        ch_target_coverages,
        ch_antitarget_coverages,
        CNVKIT_TARGET_ANTITARGET.out.target_bed,
        CNVKIT_TARGET_ANTITARGET.out.antitarget_bed,
        is_female_reference
    )

    // =========================================================================
    // Emit Results
    // =========================================================================
    emit:
    reference          = CNVKIT_REFERENCE_BUILD.out.reference
    target_bed         = CNVKIT_TARGET_ANTITARGET.out.target_bed
    antitarget_bed     = CNVKIT_TARGET_ANTITARGET.out.antitarget_bed
    alignments         = SAMTOOLS_INDEX.out.alignment
    alignment_indices  = SAMTOOLS_INDEX.out.index
    target_coverages   = CNVKIT_COVERAGE.out.target_coverage
    antitarget_coverages = CNVKIT_COVERAGE.out.antitarget_coverage
}