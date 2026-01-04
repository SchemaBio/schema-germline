/*
 * WES Trio Pipeline - Family Analysis Workflow
 *
 * 功能：全外显子测序家系（Trio）分析流程
 * 输入：先证者 + 父母的测序数据
 *
 * 流程设计：
 *   - SNV/InDel: 每个成员独立 DeepVariant → gVCF → GLnexus 联合分型
 *   - 其他模块: 以先证者为核心进行单样本分析
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
include { GLNEXUS } from '../modules/local/glnexus/main'

// Haplotype Phasing
include { WHATSHAP_PHASE } from '../modules/local/whatshap/main'

// Annotation
include { VEP } from '../modules/local/vep/main'

// STR Detection
include { EXPANSIONHUNTER } from '../modules/local/expansionhunter/main'
include { STRANGER } from '../modules/local/stranger/main'

// CNV Detection
include { CNVKIT_CALL } from '../modules/local/cnvkit/main'
include { CNVKIT_EXPORT_VCF } from '../modules/local/cnvkit/main'

// Mitochondrial Analysis
include { GATK_MUTECT2_MT } from '../modules/local/gatk/main'

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

// Genetic Model Annotation
include { GENMOD_MODELS as GENMOD_MODELS_SNP } from '../modules/local/genmod/main'
include { GENMOD_MODELS as GENMOD_MODELS_STR } from '../modules/local/genmod/main'
include { GENMOD_MODELS as GENMOD_MODELS_MT } from '../modules/local/genmod/main'
include { GENMOD_MODELS as GENMOD_MODELS_CNV } from '../modules/local/genmod/main'

// UPD Detection
include { UPD } from '../modules/local/upd/main'

// ============================================================================
// Workflow Definition
// ============================================================================

workflow WES_TRIO {

    take:
    ch_trio_reads          // Channel of [meta, [read1, read2]] - 包含 proband/father/mother
    ch_fasta
    ch_fasta_fai
    ch_bwa_amb
    ch_bwa_ann
    ch_bwa_bwt
    ch_bwa_pac
    ch_bwa_sa
    ch_bwamem2_0123
    ch_bwamem2_bwt2bit64
    use_bwamem2
    ch_target_bed
    ch_antitarget_bed
    ch_cnv_reference
    ch_variant_catalog
    ch_vep_cache
    ch_vep_plugins
    ch_vep_clinvar
    ch_vep_intervar
    ch_ped_file            // PED 家系文件 (必需)
    ch_gnotate
    ch_slivar_js
    slivar_info_expr
    slivar_sample_expr
    slivar_trio_expr       // trio 特有表达式
    slivar_family_expr
    slivar_pass_only
    ch_svdb_db
    svdb_in_occ
    svdb_in_frq
    svdb_out_occ
    svdb_out_frq
    svdb_bnd_distance
    svdb_overlap
    ch_stranger_repeats
    stranger_family_id
    ch_mt_fasta
    ch_mt_fasta_fai
    ch_mt_dict
    ch_mt_intervals
    genome_assembly

    main:
    ch_versions = Channel.empty()

    // =========================================================================
    // Step 1: FASTP - 所有成员质控
    // =========================================================================
    FASTP(ch_trio_reads)
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    // =========================================================================
    // Step 2: BWA - 所有成员比对
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
    // Step 3: MarkDuplicates - 所有成员去重
    // =========================================================================
    GATK_MARKDUPLICATES(ch_alignment, ch_fasta, 'cram')
    ch_versions = ch_versions.mix(GATK_MARKDUPLICATES.out.versions.first())

    // =========================================================================
    // Step 4: QC - 所有成员性别检查
    // =========================================================================
    SEX_CHECK(GATK_MARKDUPLICATES.out.alignment, ch_fasta, ch_fasta_fai, genome_assembly)
    ch_versions = ch_versions.mix(SEX_CHECK.out.versions.first())

    // =========================================================================
    // Step 5: DeepVariant - 所有成员生成 gVCF
    // =========================================================================
    DEEPVARIANT(
        GATK_MARKDUPLICATES.out.alignment,
        ch_fasta,
        ch_fasta_fai,
        'WES',
        ch_target_bed.ifEmpty(Channel.empty())
    )
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())

    // =========================================================================
    // Step 6: GLnexus - 联合分型
    // =========================================================================
    // 收集家系所有 gVCF
    ch_family_gvcfs = DEEPVARIANT.out.gvcf
        .map { meta, gvcf, tbi ->
            def family_id = meta.family_id ?: meta.id.split('_')[0]
            return [family_id, gvcf, tbi]
        }
        .groupTuple()
        .map { family_id, gvcfs, tbis ->
            def meta = [id: family_id, family_id: family_id]
            return [meta, gvcfs, tbis]
        }

    GLNEXUS(
        ch_family_gvcfs.map { meta, gvcfs, tbis -> meta },
        ch_family_gvcfs.map { meta, gvcfs, tbis -> gvcfs },
        ch_family_gvcfs.map { meta, gvcfs, tbis -> tbis },
        'DeepVariantWES'
    )
    ch_versions = ch_versions.mix(GLNEXUS.out.versions)

    // =========================================================================
    // Step 7: WhatsHap - 家系联合分型 (使用先证者 BAM)
    // =========================================================================
    // 获取先证者的 alignment
    ch_proband_alignment = GATK_MARKDUPLICATES.out.alignment
        .filter { meta, cram, crai -> meta.role == 'proband' }

    WHATSHAP_PHASE(
        GLNEXUS.out.vcf,
        ch_proband_alignment,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions)

    // =========================================================================
    // Step 8: VEP - 变异注释
    // =========================================================================
    VEP(
        WHATSHAP_PHASE.out.vcf,
        ch_fasta,
        ch_fasta_fai,
        ch_vep_cache.ifEmpty(file('NO_FILE')),
        ch_vep_plugins.ifEmpty(file('NO_FILE')),
        ch_vep_clinvar.ifEmpty(file('NO_FILE')),
        file('NO_FILE'),
        ch_vep_intervar.ifEmpty(file('NO_FILE')),
        file('NO_FILE'),
        file('NO_FILE'),
        file('NO_FILE'),
        file('NO_FILE'),
        file('NO_FILE'),
        file('NO_FILE'),
        '',
        '',
        genome_assembly
    )
    ch_versions = ch_versions.mix(VEP.out.versions)

    // =========================================================================
    // Step 9: GenMod - 遗传模型注释 (Trio 模式)
    // =========================================================================
    GENMOD_MODELS_SNP(
        VEP.out.vcf,
        ch_ped_file,
        true,   // split_variants
        true,   // phased
        false   // strict
    )
    ch_versions = ch_versions.mix(GENMOD_MODELS_SNP.out.versions)

    // =========================================================================
    // Step 10: Slivar - Trio 变异过滤
    // =========================================================================
    SLIVAR_EXPR(
        GENMOD_MODELS_SNP.out.vcf,
        ch_ped_file,
        ch_gnotate.ifEmpty(file('NO_FILE')),
        ch_slivar_js.ifEmpty(file('NO_FILE')),
        slivar_info_expr,
        slivar_sample_expr,
        slivar_trio_expr,
        slivar_family_expr,
        slivar_pass_only
    )
    ch_versions = ch_versions.mix(SLIVAR_EXPR.out.versions)

    // 复合杂合检测
    SLIVAR_COMPOUND_HETS(
        SLIVAR_EXPR.out.vcf,
        ch_ped_file,
        ch_gnotate.ifEmpty(file('NO_FILE'))
    )
    ch_versions = ch_versions.mix(SLIVAR_COMPOUND_HETS.out.versions)

    // TSV 导出
    SLIVAR_TSV(SLIVAR_EXPR.out.vcf, ch_ped_file)
    ch_versions = ch_versions.mix(SLIVAR_TSV.out.versions)

    // =========================================================================
    // 以下模块以先证者为核心进行分析
    // =========================================================================

    // =========================================================================
    // Step 11: STR - 先证者 STR 检测
    // =========================================================================
    if (!ch_variant_catalog.isEmpty()) {
        EXPANSIONHUNTER(
            ch_proband_alignment,
            ch_fasta,
            ch_fasta_fai,
            ch_variant_catalog
        )
        ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions)

        STRANGER(
            EXPANSIONHUNTER.out.vcf,
            ch_stranger_repeats.ifEmpty(file('NO_FILE')),
            stranger_family_id,
            false
        )
        ch_versions = ch_versions.mix(STRANGER.out.versions)

        GENMOD_MODELS_STR(
            STRANGER.out.vcf,
            ch_ped_file,
            true, false, false
        )
        ch_versions = ch_versions.mix(GENMOD_MODELS_STR.out.versions)
    }

    // =========================================================================
    // Step 12: MT - 先证者线粒体分析
    // =========================================================================
    if (!ch_mt_fasta.isEmpty()) {
        GATK_MUTECT2_MT(
            ch_proband_alignment,
            ch_mt_fasta,
            ch_mt_fasta_fai,
            ch_mt_dict,
            ch_mt_intervals
        )
        ch_versions = ch_versions.mix(GATK_MUTECT2_MT.out.versions)

        GENMOD_MODELS_MT(
            GATK_MUTECT2_MT.out.vcf,
            ch_ped_file,
            true, false, false
        )
        ch_versions = ch_versions.mix(GENMOD_MODELS_MT.out.versions)
    }

    // =========================================================================
    // Step 13: CNV - 先证者 CNV 检测
    // =========================================================================
    if (!ch_cnv_reference.isEmpty()) {
        CNVKIT_CALL(
            ch_proband_alignment,
            ch_fasta,
            ch_fasta_fai,
            ch_target_bed,
            ch_antitarget_bed,
            ch_cnv_reference
        )
        ch_versions = ch_versions.mix(CNVKIT_CALL.out.versions)

        CNVKIT_EXPORT_VCF(CNVKIT_CALL.out.call)
        ch_versions = ch_versions.mix(CNVKIT_EXPORT_VCF.out.versions)

        GENMOD_MODELS_CNV(
            CNVKIT_EXPORT_VCF.out.vcf,
            ch_ped_file,
            true, false, false
        )
        ch_versions = ch_versions.mix(GENMOD_MODELS_CNV.out.versions)

        // SVDB 注释
        SVDB_QUERY(
            CNVKIT_EXPORT_VCF.out.vcf,
            ch_svdb_db.ifEmpty(file('NO_FILE')),
            file('NO_FILE'),
            'custom',
            svdb_in_occ, svdb_in_frq,
            svdb_out_occ, svdb_out_frq,
            svdb_bnd_distance, svdb_overlap
        )
        ch_versions = ch_versions.mix(SVDB_QUERY.out.versions)
    }

    // =========================================================================
    // Step 14: UPD - 单亲二倍体检测 (Trio 特有)
    // =========================================================================
    // 获取家系成员 ID
    ch_proband_id = ch_trio_reads
        .filter { meta, reads -> meta.role == 'proband' }
        .map { meta, reads -> meta.id }

    ch_father_id = ch_trio_reads
        .filter { meta, reads -> meta.role == 'father' }
        .map { meta, reads -> meta.id }

    ch_mother_id = ch_trio_reads
        .filter { meta, reads -> meta.role == 'mother' }
        .map { meta, reads -> meta.id }

    UPD(
        GLNEXUS.out.vcf,
        ch_proband_id,
        ch_father_id,
        ch_mother_id
    )
    ch_versions = ch_versions.mix(UPD.out.versions)

    // =========================================================================
    // Emit Results
    // =========================================================================
    emit:
    // 所有成员
    fastp_json       = FASTP.out.json
    cram             = GATK_MARKDUPLICATES.out.alignment
    sex_check        = SEX_CHECK.out.result

    // 联合分型结果
    joint_vcf        = GLNEXUS.out.vcf
    phased_vcf       = WHATSHAP_PHASE.out.vcf
    vep_vcf          = VEP.out.vcf
    genmod_vcf       = GENMOD_MODELS_SNP.out.vcf

    // Slivar 过滤结果
    slivar_vcf       = SLIVAR_EXPR.out.vcf
    compound_vcf     = SLIVAR_COMPOUND_HETS.out.vcf
    slivar_tsv       = SLIVAR_TSV.out.tsv

    // 先证者特有结果
    str_vcf          = EXPANSIONHUNTER.out.vcf
    str_genmod_vcf   = GENMOD_MODELS_STR.out.vcf
    mt_vcf           = GATK_MUTECT2_MT.out.vcf
    mt_genmod_vcf    = GENMOD_MODELS_MT.out.vcf
    cnv_vcf          = CNVKIT_EXPORT_VCF.out.vcf
    cnv_genmod_vcf   = GENMOD_MODELS_CNV.out.vcf

    // UPD 结果
    upd_regions      = UPD.out.regions

    // 版本信息
    versions         = ch_versions
}
