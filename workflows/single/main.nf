/*
 * WES Single Sample Pipeline - Complete Workflow
 *
 * 功能：全外显子测序单样本分析流程
 * 包含：质控、比对、变异检测、注释、CNV、STR、线粒体等
 *
 * 输出目录结构 (简化版):
 * ${params.outdir}/
 * └── {sample_id}/
 *     ├── 01.QC/              # 质控报告 (fastp, coverage, sex_check, metrics)
 *     ├── 02.Alignment/       # 比对结果 (cram files)
 *     └── 03.Mutations/        # 突变结果 (all VCF/BCF files, 藏在子目录)
 *         ├── snv_indel/      # SNP/INDEL (deepvariant, whatshap, vep, genmod, plink2, bcftools, slivar)
 *         ├── mt/             # 线粒体 (mutect2_mt)
 *         ├── str/            # STR (expansionhunter, stranger)
 *         └── cnv/            # CNV (cnvkit, svdb)
 *
 * 存储策略:
 *   1. BWA 比对产生 sorted.cram 发布到 02.Alignment/
 *   2. GATK MarkDuplicates 完成后，marked.cram 也放入 02.Alignment/
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
include { SAMTOOLS_INDEX } from '../modules/local/samtools/main'
include { SEX_CHECK } from '../modules/local/sex_check/main'
include { XAMDST } from '../modules/local/xamdst/main'

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

// Sample Identity / Fingerprint
include { SAMPLE_FINGERPRINT } from '../modules/local/sample_identity/main'

// Variant Filtering
include { SLIVAR_EXPR } from '../modules/local/slivar/main'
include { SLIVAR_COMPOUND_HETS } from '../modules/local/slivar/main'
include { SLIVAR_TSV } from '../modules/local/slivar/main'

// CNV/SV Annotation
include { SVDB_QUERY } from '../modules/local/svdb/main'

// STR Annotation
include { STRANGER } from '../modules/local/stranger/main'

// Mitochondrial Analysis
include { GATK_MUTECT2_MT } from '../modules/local/gatk/main'

// Genetic Model Annotation
include { GENMOD_MODELS } from '../modules/local/genmod/main'
include { GENMOD_ANNOTATE } from '../modules/local/genmod/main'
include { GENMOD_SORT } from '../modules/local/genmod/main'
include { GENMOD_MODELS as GENMOD_MODELS_SNP } from '../modules/local/genmod/main'
include { GENMOD_MODELS as GENMOD_MODELS_STR } from '../modules/local/genmod/main'
include { GENMOD_MODELS as GENMOD_MODELS_MT } from '../modules/local/genmod/main'
include { GENMOD_MODELS as GENMOD_MODELS_CNV } from '../modules/local/genmod/main'


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

    main:
    ch_versions = Channel.empty()

    // =========================================================================
    // 默认可选参数 (在 main 块中定义)
    // =========================================================================
    ch_target_bed      = Channel.empty()      // 目标区域 BED
    ch_antitarget_bed  = Channel.empty()      // 反目标区域 BED
    ch_cnv_reference   = Channel.empty()      // CNV 参考基线
    ch_variant_catalog = Channel.empty()      // STR 位点目录
    // VEP 参数
    ch_vep_cache       = Channel.empty()      // VEP 缓存
    ch_vep_plugins     = Channel.empty()      // VEP 插件
    ch_vep_clinvar     = Channel.empty()      // ClinVar VCF
    ch_vep_intervar    = Channel.empty()      // InterVar VCF
    ch_vep_gnomad      = Channel.empty()      // gnomAD VCF
    ch_vep_dbsnp       = Channel.empty()      // dbSNP VCF
    ch_vep_alphamissense = Channel.empty()    // AlphaMissense DB
    ch_vep_evoscore2   = Channel.empty()      // EVOScore2 DB
    ch_vep_pangolin    = Channel.empty()      // Pangolin DB
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
    enable_mt_analysis = true                 // 是否开启线粒体分析（使用主参考基因组）
    // Sample Fingerprint 参数
    ch_snp_panel = Channel.value(file('NO_FILE'))  // 样本识别 SNP 位点 VCF
    // GenMod 参数
    ch_genmod_ped = Channel.value(file('NO_FILE'))  // PED 家系文件
    genmod_split_variants = true              // 拆分多等位基因
    genmod_phased = true                      // 输入已定相
    genmod_strict = false                     // 严格模式
    genome_assembly    = 'GRCh38'              // 基因组版本

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
    // Step 3b: SAMTOOLS_INDEX - Create BAM/CRAM Index
    // =========================================================================
    SAMTOOLS_INDEX(GATK_MARKDUPLICATES.out.alignment)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // =========================================================================
    // Step 4: SEX_CHECK - Sex Verification (based on SRY gene coverage)
    // =========================================================================
    SEX_CHECK(SAMTOOLS_INDEX.out.alignment, genome_assembly)
    ch_versions = ch_versions.mix(SEX_CHECK.out.versions.first())

    // =========================================================================
    // Step 5: XAMDST - Coverage Statistics
    // =========================================================================
    if (!ch_target_bed.isEmpty()) {
        XAMDST(SAMTOOLS_INDEX.out.alignment, ch_target_bed, Channel.empty())
        ch_versions = ch_versions.mix(XAMDST.out.versions.first())
    } else {
        ch_xamdst_out = Channel.empty()
    }

    // =========================================================================
    // Step 6: DEEPVARIANT - SNP/INDEL Variant Calling
    // =========================================================================
    DEEPVARIANT(
        SAMTOOLS_INDEX.out.alignment,
        ch_fasta,
        ch_fasta_fai,
        'WES',
        !ch_target_bed.isEmpty() ? ch_target_bed : Channel.empty()
    )
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())

    // =========================================================================
    // Step 6b: SAMPLE_FINGERPRINT - Sample Identity Verification
    // =========================================================================
    if (!ch_snp_panel.isEmpty()) {
        SAMPLE_FINGERPRINT(
            DEEPVARIANT.out.vcf,
            ch_snp_panel
        )
        ch_versions = ch_versions.mix(SAMPLE_FINGERPRINT.out.versions.first())
    } else {
        ch_fingerprint_vcf = Channel.empty()
        ch_fingerprint_txt = Channel.empty()
    }

    // =========================================================================
    // Step 6c: GATK MUTECT2 - Mitochondrial Variant Calling
    // =========================================================================
    if (enable_mt_analysis) {
        GATK_MUTECT2_MT(
            SAMTOOLS_INDEX.out.alignment,
            ch_fasta,
            ch_fasta_fai
        )
        ch_versions = ch_versions.mix(GATK_MUTECT2_MT.out.versions.first())

        // =========================================================================
        // Step 6c: GENMOD - MT Genetic Model Annotation
        // =========================================================================
        GENMOD_MODELS_MT(
            GATK_MUTECT2_MT.out.vcf,
            ch_genmod_ped,
            genmod_split_variants,
            false,  // MT 变异通常未定相
            genmod_strict
        )
        ch_versions = ch_versions.mix(GENMOD_MODELS_MT.out.versions.first())
    } else {
        ch_mt_vcf = Channel.empty()
    }

    // =========================================================================
    // Step 7: WHATSHAP - Haplotype Phasing
    // =========================================================================
    WHATSHAP_PHASE(
        DEEPVARIANT.out.vcf,
        SAMTOOLS_INDEX.out.alignment,
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
        // Custom VCF databases
        !ch_vep_clinvar.isEmpty() ? ch_vep_clinvar : Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        !ch_vep_intervar.isEmpty() ? ch_vep_intervar : Channel.value(file('NO_FILE')),
        !ch_vep_gnomad.isEmpty() ? ch_vep_gnomad : Channel.value(file('NO_FILE')),
        !ch_vep_dbsnp.isEmpty() ? ch_vep_dbsnp : Channel.value(file('NO_FILE')),
        // Plugin databases
        !ch_vep_alphamissense.isEmpty() ? ch_vep_alphamissense : Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        !ch_vep_evoscore2.isEmpty() ? ch_vep_evoscore2 : Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        !ch_vep_pangolin.isEmpty() ? ch_vep_pangolin : Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        // Other databases
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        Channel.value(file('NO_FILE')),
        '',
        '',
        genome_assembly
    )
    ch_versions = ch_versions.mix(VEP.out.versions.first())

    // =========================================================================
    // Step 8b: GENMOD - Genetic Model Annotation (在 VEP 之后)
    // =========================================================================
    GENMOD_MODELS_SNP(
        VEP.out.vcf,
        ch_genmod_ped,
        genmod_split_variants,
        genmod_phased,
        genmod_strict
    )
    ch_versions = ch_versions.mix(GENMOD_MODELS_SNP.out.versions.first())

    // =========================================================================
    // Step 9: EXPANSIONHUNTER - STR Detection
    // =========================================================================
    if (!ch_variant_catalog.isEmpty()) {
        EXPANSIONHUNTER(
            SAMTOOLS_INDEX.out.alignment,
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

        // =========================================================================
        // Step 9c: GENMOD - STR Genetic Model Annotation (在 Stranger 之后)
        // =========================================================================
        GENMOD_MODELS_STR(
            STRANGER.out.vcf,
            ch_genmod_ped,
            genmod_split_variants,
            genmod_phased,
            genmod_strict
        )
        ch_versions = ch_versions.mix(GENMOD_MODELS_STR.out.versions.first())
    } else {
        ch_expansionhunter_out = Channel.empty()
    }

    // =========================================================================
    // Step 10: CNVKIT - CNV Detection
    // =========================================================================
    if (!ch_cnv_reference.isEmpty() && !ch_target_bed.isEmpty() && !ch_antitarget_bed.isEmpty()) {
        CNVKIT_CALL(
            SAMTOOLS_INDEX.out.alignment,
            ch_fasta,
            ch_fasta_fai,
            ch_target_bed,
            ch_antitarget_bed,
            ch_cnv_reference
        )
        ch_versions = ch_versions.mix(CNVKIT_CALL.out.versions.first())

        // Export CNV to VCF
        CNVKIT_EXPORT_VCF(CNVKIT_CALL.out.call)
        ch_versions = ch_versions.mix(CNVKIT_EXPORT_VCF.out.versions.first())

        // =========================================================================
        // Step 10b: GENMOD - CNV Genetic Model Annotation
        // =========================================================================
        GENMOD_MODELS_CNV(
            CNVKIT_EXPORT_VCF.out.vcf,
            ch_genmod_ped,
            genmod_split_variants,
            false,  // CNV 变异未定相
            genmod_strict
        )
        ch_versions = ch_versions.mix(GENMOD_MODELS_CNV.out.versions.first())
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
    // Step 13: SLIVAR - Variant Filtering (使用 VEP 注释后的 VCF)
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
    cram = SAMTOOLS_INDEX.out.alignment
    markdup_metrics = GATK_MARKDUPLICATES.out.metrics

    // QC
    sex_check = SEX_CHECK.out.result
    inferred_sex = SEX_CHECK.out.sex
    coverage = !ch_xamdst_out.isEmpty() ? ch_xamdst_out.coverage_stat : Channel.empty()
    coverage_stat = !ch_xamdst_out.isEmpty() ? ch_xamdst_out.coverage_stat : Channel.empty()
    coverage_json = !ch_xamdst_out.isEmpty() ? ch_xamdst_out.json : Channel.empty()
    chromosome_stat = !ch_xamdst_out.isEmpty() ? ch_xamdst_out.chromosome_stat : Channel.empty()
    region_stat = !ch_xamdst_out.isEmpty() ? ch_xamdst_out.region_stat : Channel.empty()

    // Variants
    deepvariant_vcf = DEEPVARIANT.out.vcf
    deepvariant_gvcf = DEEPVARIANT.out.gvcf
    deepvariant_report = DEEPVARIANT.out.report

    // Sample Fingerprint
    fingerprint_vcf = !ch_snp_panel.isEmpty() ? ch_fingerprint_vcf : Channel.empty()
    fingerprint_txt = !ch_snp_panel.isEmpty() ? ch_fingerprint_txt : Channel.empty()

    // Mitochondrial variants
    mt_vcf = GATK_MUTECT2_MT.out.vcf
    mt_stats = GATK_MUTECT2_MT.out.stats

    // MT GenMod annotated
    mt_genmod_vcf = GENMOD_MODELS_MT.out.vcf

    // Phased variants (WhatsHap)
    whatshap_vcf = WHATSHAP_PHASE.out.vcf

    // Annotated variants
    vep_vcf = VEP.out.vcf

    // GenMod annotated variants (SNP/INDEL)
    genmod_vcf = GENMOD_MODELS_SNP.out.vcf

    // STR
    str_vcf = !ch_expansionhunter_out.isEmpty() ? ch_expansionhunter_out.vcf : Channel.empty()
    str_json = !ch_expansionhunter_out.isEmpty() ? ch_expansionhunter_out.json : Channel.empty()

    // STR annotated (stranger)
    str_annotated_vcf = STRANGER.out.vcf

    // STR GenMod annotated
    str_genmod_vcf = !ch_expansionhunter_out.isEmpty() ? GENMOD_MODELS_STR.out.vcf : Channel.empty()

    // CNV
    cnv_cnr = !ch_cnvkit_out.isEmpty() ? ch_cnvkit_out.cnr : Channel.empty()
    cnv_cns = !ch_cnvkit_out.isEmpty() ? ch_cnvkit_out.cns : Channel.empty()
    cnv_call = !ch_cnvkit_out.isEmpty() ? ch_cnvkit_out.call : Channel.empty()
    cnv_vcf = !ch_cnvkit_vcf.isEmpty() ? ch_cnvkit_vcf.vcf : Channel.empty()

    // CNV GenMod annotated
    cnv_genmod_vcf = GENMOD_MODELS_CNV.out.vcf

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
