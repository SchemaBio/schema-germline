/*
 * WES Single Sample Pipeline - Complete Workflow
 *
 * 功能：全外显子测序单样本分析流程
 * 包含：质控、比对、变异检测、注释、CNV、STR 等
 *
 * 输出目录结构 (每个样本一个目录):
 * results/
 * └── {sample_id}/
 *     ├── 01_fastp/          # 质控报告
 *     ├── 02_alignment/      # 比对结果 (sorted.cram -> marked.cram)
 *     ├── 03_markdup/        # 去重统计
 *     ├── 04_coverage/       # 覆盖度统计
 *     ├── 05_sex_check/      # 性别检查
 *     ├── 06_deepvariant/    # SNP/INDEL 变异检测
 *     ├── 07_whatshap/       # 单倍型分相 (WhatsHap)
 *     ├── 07_vep/            # 变异注释
 *     ├── 08_str/            # STR 检测 (ExpansionHunter)
 *     ├── 09_str/            # STR 注释 (stranger)
 *     ├── 10_cnv/            # CNV 检测
 *     ├── 11_plink2/         # 基因型分析
 *     ├── 12_roh/            # ROH 分析
 *     ├── 13_phasing/        # Phasing 统计
 *     ├── 14_slivar/         # 变异过滤 (slivar)
 *     └── 15_svdb/           # CNV/SV 注释 (SVDB)
 *
 * 存储策略:
 *   1. BWA 比对产生 sorted.cram 发布到 02_alignment/
 *   2. GATK MarkDuplicates 完成后，marked.cram 也放入 02_alignment/
 *   3. 然后自动删除 sorted.cram，节约空间
 */

// ============================================================================
// Include Modules
// ============================================================================

// Quality Control
include { FASTP } from '../modules/local/fastp/main'

// Alignment
include { BWA_MEM } from '../modules/local/bwa_mem2/main'
include { BWA_MEM2 } from '../modules/local/bwa_mem2/main'

// Post-Alignment Processing
include { GATK_MARKDUPLICATES } from '../modules/local/gatk/main'
include { SEX_CHECK } from '../modules/local/sex_check/main'
include { BAMDST } from '../modules/local/bamdst/main'

// Variant Calling
include { DEEPVARIANT } from '../modules/local/deepvariant/main'

// Haplotype Phasing
include { WHATSHAP_PHASE } from '../modules/local/whatshap/main'

// Annotation
include { VEP } from '../modules/local/vep/main'

// STR Detection
include { EXPANSIONHUNTER } from '../modules/local/expansionhunter/main'

// CNV Detection
include { CNVKIT_CALL } from '../modules/local/cnvkit/main'
include { CNVKIT_EXPORT_VCF } from '../modules/local/cnvkit/main'

// Genetic Analysis
include { PLINK2_VCF_TO_PLINK } from '../modules/local/plink2/main'
include { BCFTOOLS_ROH } from '../modules/local/bcftools/main'
include { BCFTOOLS_STATS } from '../modules/local/bcftools/main'

// Variant Filtering
include { SLIVAR_EXPR } from '../modules/local/slivar/main'
include { SLIVAR_COMPOUND_HETS } from '../modules/local/slivar/main'
include { SLIVAR_TSV } from '../modules/local/slivar/main'

// CNV/SV Annotation
include { SVDB_QUERY } from '../modules/local/svdb/main'

// STR Annotation
include { STRANGER } from '../modules/local/stranger/main'


// ============================================================================
// Workflow Definition
// ============================================================================

workflow WES_SINGLE {

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
    // 可选参数
    ch_target_bed      = Channel.empty()      // 目标区域 BED
    ch_antitarget_bed  = Channel.empty()      // 反目标区域 BED
    ch_cnv_reference   = Channel.empty()      // CNV 参考基线
    ch_variant_catalog = Channel.empty()      // STR 位点目录
    ch_vep_cache       = Channel.empty()      // VEP 缓存
    ch_vep_plugins     = Channel.empty()      // VEP 插件
    ch_vep_clinvar     = Channel.empty()      // ClinVar VCF
    ch_vep_intervar    = Channel.empty()      // InterVar VCF
    ch_ped_file        = Channel.value(file('NO_FILE'))  // PED 系谱文件
    ch_gnotate         = Channel.value(file('NO_FILE'))  // gnomad gnotate 文件
    ch_slivar_js       = Channel.value(file('NO_FILE'))  // slivar-functions.js
    slivar_info_expr   = ''                    // slivar INFO 过滤表达式
    slivar_sample_expr = ''                    // slivar 样本表达式
    slivar_trio_expr   = ''                    // slivar 三联体表达式
    slivar_family_expr = ''                    // slivar 家族表达式
    slivar_pass_only   = false                 // slivar 只输出 PASS 变体
    // SVDB 参数
    ch_svdb_db         = Channel.value(file('NO_FILE'))  // SVDB 数据库 VCF
    svdb_in_occ        = ''                    // 数据库计数标签
    svdb_in_frq        = ''                    // 数据库频率标签
    svdb_out_occ       = ''                    // 输出计数标签名
    svdb_out_frq       = ''                    // 输出频率标签名
    svdb_bnd_distance  = 2500                  // 断点距离阈值
    svdb_overlap       = 0.8                   // 重叠比例阈值
    // Stranger 参数
    ch_stranger_repeats = Channel.value(file('NO_FILE'))  // STR 重复定义 JSON
    stranger_family_id = ''                   // 家族 ID
    genome_assembly    = 'GRCh38'              // 基因组版本

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
    // Step 3: GATK MarkDuplicates - Mark PCR Duplicates
    // =========================================================================
    GATK_MARKDUPLICATES(ch_alignment, ch_fasta, 'cram')
    ch_versions = ch_versions.mix(GATK_MARKDUPLICATES.out.versions.first())

    // =========================================================================
    // Step 4: SEX_CHECK - Sex Verification (based on SRY gene coverage)
    // =========================================================================
    SEX_CHECK(GATK_MARKDUPLICATES.out.alignment, ch_fasta, ch_fasta_fai, genome_assembly)
    ch_versions = ch_versions.mix(SEX_CHECK.out.versions.first())

    // =========================================================================
    // Step 5: BAMDST - Coverage Statistics
    // =========================================================================
    if (!ch_target_bed.isEmpty()) {
        BAMDST(GATK_MARKDUPLICATES.out.alignment, ch_target_bed, Channel.empty())
        ch_versions = ch_versions.mix(BAMDST.out.versions.first())
    } else {
        ch_bamdst_out = Channel.empty()
    }

    // =========================================================================
    // Step 6: DEEPVARIANT - SNP/INDEL Variant Calling
    // =========================================================================
    DEEPVARIANT(
        GATK_MARKDUPLICATES.out.alignment,
        ch_fasta,
        ch_fasta_fai,
        'WES',
        !ch_target_bed.isEmpty() ? ch_target_bed : Channel.empty()
    )
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())

    // =========================================================================
    // Step 7: WHATSHAP - Haplotype Phasing
    // =========================================================================
    WHATSHAP_PHASE(
        DEEPVARIANT.out.vcf,
        GATK_MARKDUPLICATES.out.alignment,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions.first())

    // =========================================================================
    // Step 8: VEP - Variant Annotation
    // =========================================================================
    VEP(
        WHATSHAP_PHASE.out.vcf,
        ch_fasta,
        ch_fasta_fai,
        !ch_vep_cache.isEmpty() ? ch_vep_cache : Channel.value(file('NO_FILE')),
        !ch_vep_plugins.isEmpty() ? ch_vep_plugins : Channel.value(file('NO_FILE')),
        !ch_vep_clinvar.isEmpty() ? ch_vep_clinvar : Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        !ch_vep_intervar.isEmpty() ? ch_vep_intervar : Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        '',
        '',
        genome_assembly
    )
    ch_versions = ch_versions.mix(VEP.out.versions.first())

    // =========================================================================
    // Step 9: EXPANSIONHUNTER - STR Detection
    // =========================================================================
    if (!ch_variant_catalog.isEmpty()) {
        EXPANSIONHUNTER(
            GATK_MARKDUPLICATES.out.alignment,
            ch_fasta,
            ch_fasta_fai,
            ch_variant_catalog
        )
        ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions.first())

        // =========================================================================
        // Step 9b: STRANGER - STR Annotation (在 ExpansionHunter 之后)
        // =========================================================================
        STRANGER(
            EXPANSIONHUNTER.out.vcf,
            ch_stranger_repeats,
            stranger_family_id,
            false
        )
        ch_versions = ch_versions.mix(STRANGER.out.versions.first())
    } else {
        ch_expansionhunter_out = Channel.empty()
    }

    // =========================================================================
    // Step 10: CNVKIT - CNV Detection
    // =========================================================================
    if (!ch_cnv_reference.isEmpty() && !ch_target_bed.isEmpty() && !ch_antitarget_bed.isEmpty()) {
        CNVKIT_CALL(
            GATK_MARKDUPLICATES.out.alignment,
            ch_fasta,
            ch_fasta_fai,
            ch_target_bed,
            ch_antitarget_bed,
            ch_cnv_reference
        )
        ch_versions = ch_versions.mix(CNVKIT_CALL.out.versions.first())

        // Export CNV to VCF
        CNVKIT_EXPORT_VCF(CNVKIT_CALL.out.call)
    } else {
        ch_cnvkit_out = Channel.empty()
        ch_cnvkit_vcf = Channel.empty()
    }

    // =========================================================================
    // Step 11: PLINK2 - Genotype Data Analysis (使用 Phased VCF)
    // =========================================================================
    PLINK2_VCF_TO_PLINK(WHATSHAP_PHASE.out.vcf, ch_fasta_fai)
    ch_versions = ch_versions.mix(PLINK2_VCF_TO_PLINK.out.versions.first())

    // =========================================================================
    // Step 12: BCFTOOLS - ROH & VCF Stats (使用 Phased VCF)
    // =========================================================================
    BCFTOOLS_ROH(WHATSHAP_PHASE.out.vcf, ch_fasta, ch_fasta_fai)
    ch_versions = ch_versions.mix(BCFTOOLS_ROH.out.versions.first())

    BCFTOOLS_STATS(WHATSHAP_PHASE.out.vcf, ch_fasta, ch_fasta_fai)
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())

    // =========================================================================
    // Step 13: VEP - Variant Annotation
    // =========================================================================
    VEP(
        WHATSHAP_PHASE.out.vcf,
        ch_fasta,
        ch_fasta_fai,
        !ch_vep_cache.isEmpty() ? ch_vep_cache : Channel.value(file('NO_FILE')),
        !ch_vep_plugins.isEmpty() ? ch_vep_plugins : Channel.value(file('NO_FILE')),
        !ch_vep_clinvar.isEmpty() ? ch_vep_clinvar : Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        !ch_vep_intervar.isEmpty() ? ch_vep_intervar : Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        '',
        '',
        genome_assembly
    )
    ch_versions = ch_versions.mix(VEP.out.versions.first())

    // =========================================================================
    // Step 14: SLIVAR - Variant Filtering (使用 VEP 注释后的 VCF)
    // =========================================================================
    SLIVAR_EXPR(
        VEP.out.vcf,
        ch_ped_file,
        ch_gnotate,
        ch_slivar_js,
        slivar_info_expr,
        slivar_sample_expr,
        slivar_trio_expr,
        slivar_family_expr,
        slivar_pass_only
    )
    ch_versions = ch_versions.mix(SLIVAR_EXPR.out.versions.first())

    // SLIVAR compound-hets detection (可选)
    SLIVAR_COMPOUND_HETS(
        SLIVAR_EXPR.out.vcf,
        ch_ped_file,
        ch_gnotate
    )
    ch_versions = ch_versions.mix(SLIVAR_COMPOUND_HETS.out.versions.first())

    // SLIVAR TSV export (可选)
    SLIVAR_TSV(
        SLIVAR_EXPR.out.vcf,
        ch_ped_file
    )
    ch_versions = ch_versions.mix(SLIVAR_TSV.out.versions.first())

    // =========================================================================
    // Step 15: SVDB - CNV/SV Annotation (使用 CNV VCF)
    // =========================================================================
    SVDB_QUERY(
        !ch_cnvkit_vcf.isEmpty() ? ch_cnvkit_vcf.vcf : Channel.empty(),
        ch_svdb_db,
        Channel.value(file('NO_FILE')),
        'custom',
        svdb_in_occ,
        svdb_in_frq,
        svdb_out_occ,
        svdb_out_frq,
        svdb_bnd_distance,
        svdb_overlap
    )
    ch_versions = ch_versions.mix(SVDB_QUERY.out.versions.first())

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

    // Alignment
    cram = GATK_MARKDUPLICATES.out.alignment
    cram_index = GATK_MARKDUPLICATES.out.alignment.map { it[2] }
    markdup_metrics = GATK_MARKDUPLICATES.out.metrics

    // QC
    sex_check = SEX_CHECK.out.result
    inferred_sex = SEX_CHECK.out.sex
    coverage = !ch_bamdst_out.isEmpty() ? ch_bamdst_out.coverage_stat : Channel.empty()
    coverage_stat = !ch_bamdst_out.isEmpty() ? ch_bamdst_out.coverage_stat : Channel.empty()
    chromosome_stat = !ch_bamdst_out.isEmpty() ? ch_bamdst_out.chromosome_stat : Channel.empty()
    region_stat = !ch_bamdst_out.isEmpty() ? ch_bamdst_out.region_stat : Channel.empty()

    // Variants
    deepvariant_vcf = DEEPVARIANT.out.vcf
    deepvariant_gvcf = DEEPVARIANT.out.gvcf
    deepvariant_report = DEEPVARIANT.out.report

    // Phased variants (WhatsHap)
    whatshap_vcf = WHATSHAP_PHASE.out.vcf

    // Annotated variants
    vep_vcf = VEP.out.vcf

    // STR
    str_vcf = !ch_expansionhunter_out.isEmpty() ? ch_expansionhunter_out.vcf : Channel.empty()
    str_json = !ch_expansionhunter_out.isEmpty() ? ch_expansionhunter_out.json : Channel.empty()

    // STR annotated (stranger)
    str_annotated_vcf = STRANGER.out.vcf

    // CNV
    cnv_cnr = !ch_cnvkit_out.isEmpty() ? ch_cnvkit_out.cnr : Channel.empty()
    cnv_cns = !ch_cnvkit_out.isEmpty() ? ch_cnvkit_out.cns : Channel.empty()
    cnv_call = !ch_cnvkit_out.isEmpty() ? ch_cnvkit_out.call : Channel.empty()
    cnv_vcf = !ch_cnvkit_vcf.isEmpty() ? ch_cnvkit_vcf.vcf : Channel.empty()

    // Plink2
    plink2_pgen = PLINK2_VCF_TO_PLINK.out.plink2.map { it[0] }
    plink2_pvar = PLINK2_VCF_TO_PLINK.out.plink2.map { it[1] }
    plink2_psam = PLINK2_VCF_TO_PLINK.out.plink2.map { it[2] }

    // ROH & Stats
    roh = BCFTOOLS_ROH.out.roh
    vcf_stats = BCFTOOLS_STATS.out.stats

    // Slivar filtered variants
    slivar_vcf = SLIVAR_EXPR.out.vcf
    slivar_compound_vcf = SLIVAR_COMPOUND_HETS.out.vcf
    slivar_tsv = SLIVAR_TSV.out.tsv

    // SVDB annotated CNV
    svdb_vcf = SVDB_QUERY.out.vcf

    // Versions
    versions = ch_versions
}
