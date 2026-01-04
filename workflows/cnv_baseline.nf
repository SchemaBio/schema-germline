/*
 * CNV Baseline Pipeline
 *
 * 功能：使用多个正常样本构建 CNVkit 参考基线
 * 输入：多个正常样本的 BAM/CRAM 文件
 * 输出：reference.cnn 基线文件
 *
 * 使用方法：
 *   nextflow run main.nf -entry CNV_BASELINE \
 *       --input samples.csv \
 *       --fasta /path/to/reference.fa \
 *       --target_bed /path/to/targets.bed \
 *       -profile docker
 *
 * 输入 CSV 格式：
 *   sample_id,bam,bai
 *   normal1,/path/to/normal1.cram,/path/to/normal1.cram.crai
 *   normal2,/path/to/normal2.bam,/path/to/normal2.bam.bai
 *
 * 说明：
 *   - 至少需要 1 个正常样本，推荐 >= 5 个以获得稳定基线
 *   - 支持 BAM 和 CRAM 格式
 */

// ============================================================================
// Include Modules
// ============================================================================

include { CNVKIT_TARGET } from '../modules/local/cnvkit/main'
include { CNVKIT_ANTITARGET } from '../modules/local/cnvkit/main'
include { CNVKIT_COVERAGE } from '../modules/local/cnvkit/main'

// ============================================================================
// Workflow Definition
// ============================================================================

workflow CNV_BASELINE {

    take:
    ch_alignments      // Channel of [meta, bam, bai] tuples
    ch_fasta           // 参考基因组
    ch_fasta_fai       // 参考基因组索引
    ch_target_bed      // 目标区域 BED (捕获区域)
    ch_annotate        // 注释文件 refFlat.txt (可选)
    ch_access_bed      // 可访问区域 BED (可选)

    main:
    ch_versions = Channel.empty()

    // =========================================================================
    // Step 1: 准备 Target BED
    // =========================================================================
    CNVKIT_TARGET(
        ch_target_bed,
        ch_annotate.ifEmpty(file('NO_FILE'))
    )
    ch_versions = ch_versions.mix(CNVKIT_TARGET.out.versions)

    // =========================================================================
    // Step 2: 生成 Antitarget BED
    // =========================================================================
    CNVKIT_ANTITARGET(
        CNVKIT_TARGET.out.target,
        ch_access_bed.ifEmpty(file('NO_FILE'))
    )
    ch_versions = ch_versions.mix(CNVKIT_ANTITARGET.out.versions)

    // =========================================================================
    // Step 3: 计算每个样本的覆盖度
    // =========================================================================
    CNVKIT_COVERAGE(
        ch_alignments,
        ch_fasta,
        ch_fasta_fai,
        CNVKIT_TARGET.out.target,
        CNVKIT_ANTITARGET.out.antitarget
    )
    ch_versions = ch_versions.mix(CNVKIT_COVERAGE.out.versions.first())

    // =========================================================================
    // Step 4: 收集覆盖度文件并构建基线
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
    target_bed     = CNVKIT_TARGET.out.target
    antitarget_bed = CNVKIT_ANTITARGET.out.antitarget
    versions       = ch_versions
}

// ============================================================================
// 构建基线进程
// ============================================================================

process CNVKIT_REFERENCE_BUILD {
    tag "reference"
    label 'process_medium'

    input:
    path target_cnns       // *.targetcoverage.cnn 文件
    path antitarget_cnns   // *.antitargetcoverage.cnn 文件
    path fasta

    output:
    path "reference.cnn", emit: reference
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    cnvkit.py reference \\
        ${target_cnns} \\
        ${antitarget_cnns} \\
        --fasta ${fasta} \\
        --output reference.cnn \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed 's/cnvkit //')
    END_VERSIONS
    """
}
