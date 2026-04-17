/*
 * CNV Baseline Pipeline
 *
 * 功能：使用多个正常样本构建 CNVkit 参考基线
 * 输入：多个正常样本的 BAM 文件 + 预制的 target/antitarget BED
 * 输出：reference.cnn 基线文件
 *
 * 使用方法：
 *   nextflow run main.nf -entry CNV_BASELINE \
 *       --input samples.csv \
 *       --fasta /path/to/reference.fa \
 *       --target_bed /path/to/capture.bed \
 *       -profile docker
 *
 * 输入 CSV 格式：
 *   sample_id,bam,bai
 *   normal1,/path/to/normal1.bam,/path/to/normal1.bam.bai
 *   normal2,/path/to/normal2.bam,/path/to/normal2.bam.bai
 *
 * 说明：
 *   - 至少需要 1 个正常样本，推荐 >= 5 个以获得稳定基线
 *   - 支持 BAM 格式
 *   - 使用预制的带注释 target/antitarget BED，与用户 BED 取交集
 */

// ============================================================================
// Include Modules
// ============================================================================

include { CNVKIT_INTERSECT } from '../modules/local/cnvkit/main'
include { CNVKIT_COVERAGE } from '../modules/local/cnvkit/main'
include { CNVKIT_REFERENCE_BUILD } from '../modules/local/cnvkit/main'

// ============================================================================
// Workflow Definition
// ============================================================================

workflow CNV_BASELINE {

    take:
    ch_alignments      // Channel of [meta, bam, bai] tuples
    ch_fasta           // 参考基因组
    ch_fasta_fai       // 参考基因组索引
    ch_target_bed      // 用户捕获区域 BED
    ch_prebuilt_target     // 预制 target BED (含注释)
    ch_prebuilt_antitarget // 预制 antitarget BED (含注释)

    main:
    ch_versions = Channel.empty()

    // =========================================================================
    // Step 1: 用户 BED 与预制 target/antitarget 取交集
    // =========================================================================
    CNVKIT_INTERSECT(
        ch_target_bed,
        ch_prebuilt_target,
        ch_prebuilt_antitarget
    )
    ch_versions = ch_versions.mix(CNVKIT_INTERSECT.out.versions)

    // =========================================================================
    // Step 2: 计算每个样本的覆盖度
    // =========================================================================
    CNVKIT_COVERAGE(
        ch_alignments,
        ch_fasta,
        ch_fasta_fai,
        CNVKIT_INTERSECT.out.target,
        CNVKIT_INTERSECT.out.antitarget
    )
    ch_versions = ch_versions.mix(CNVKIT_COVERAGE.out.versions.first())

    // =========================================================================
    // Step 3: 收集覆盖度文件并构建基线
    // =========================================================================
    ch_target_cnns = CNVKIT_COVERAGE.out.target_coverage
        .map { meta, cnn -> cnn }
        .collect()

    ch_antitarget_cnns = CNVKIT_COVERAGE.out.antitarget_coverage
        .map { meta, cnn -> cnn }
        .collect()

    CNVKIT_REFERENCE_BUILD(
        ch_target_cnns,
        ch_antitarget_cnns,
        ch_fasta
    )
    ch_versions = ch_versions.mix(CNVKIT_REFERENCE_BUILD.out.versions)

    // =========================================================================
    // Emit Results
    // =========================================================================
    emit:
    reference      = CNVKIT_REFERENCE_BUILD.out.reference
    target_bed     = CNVKIT_INTERSECT.out.target
    antitarget_bed = CNVKIT_INTERSECT.out.antitarget
    versions       = ch_versions
}
