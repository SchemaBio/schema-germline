/*
 * WES Trio Pipeline
 *
 * 功能：家系全外显子测序 (WES) 分析流程
 * 输入：家系成员 FASTQ 文件 (父母子 Trio) + 参考基因组 + 捕获区域 BED + 配置参数
 * 输出：家系联合变异检测结果 + 先证者专属分析 (CNV/STR/MEI) + 亲缘验证 + ROH/UPD
 *
 * 流程设计:
 *   1. 所有成员进行 QC + 比对 + DeepVariant (gVCF)
 *   2. GLnexus 合并家系 gVCF 为 joint VCF
 *   3. Peddy 验证亲缘关系
 *   4. 对合并 VCF 进行定相 + VEP 注释
 *   5. 先证者专属分析: CNV, STR, MEI, 线粒体
 *   6. Trio ROH/UPD 检测 (父母子三人)
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
include { CNVKIT_BATCH } from '../../modules/cnvkit/main'
include { EXPANSIONHUNTER } from '../../modules/expansionhunter/main'
include { TIEA_WES } from '../../modules/tiea-wes/main'
include { WHATSHAP_PHASE } from '../../modules/whatshap/main'
include { VEP_ANNOTATE; VEP_MT; VEP_MEI } from '../../modules/vep/main'
include { STRANGER_ANNOTATE; STRANGER_FILTER } from '../../modules/stranger/main'
include { AUTOMAP_ROH; AUTOMAP_UPD } from '../../modules/automap/main'
include { GLNEXUS_JOINTCALL_TRIO } from '../../modules/glnexus/main'
include { PEDDY_FAMILY } from '../../modules/peddy/main'

// ============================================================================
// Workflow Definition
// ============================================================================

workflow WES_TRIO_PIPELINE {

    take:
    ch_family_reads           // Channel of [family_id, member_role, sample_id, [read1, read2]]
    ch_ped_file               // PED 文件
    ch_fasta                  // 参考基因组 FASTA
    ch_fasta_fai              // 参考基因组索引
    ch_fasta_dict             // 参考基因组字典
    ch_target_bed             // 捕获区域 BED
    ch_mt_bed                 // 线粒体区域 BED
    ch_snp_positions          // SNP 位点文件
    ch_cnvkit_reference       // CNVkit 基线
    ch_cnvkit_access_bed      // CNVkit access BED
    ch_expansionhunter_catalog // ExpansionHunter catalog
    ch_vep_db_files           // VEP 数据库文件 (可选)
    val family_config         // 配置参数 Map

    main:
    // =========================================================================
    // 解析配置参数
    // =========================================================================
    def use_gpu = family_config.use_gpu ?: false
    def genome_assembly = family_config.genome_assembly ?: 'GRCh38'
    def sex_check_threshold = family_config.sex_check_threshold ?: 10
    def baf_max_depth = family_config.baf_max_depth ?: 200
    def deepvariant_model = family_config.deepvariant_model ?: 'WES'
    def deepvariant_shards = family_config.deepvariant_shards
    def glnexus_config = family_config.glnexus_config ?: 'deepvariant'
    def cnvkit_annotate = family_config.cnvkit_annotate ?: true
    def cnvkit_split_size = family_config.cnvkit_split_size ?: 5000
    def cnvkit_min_target_size = family_config.cnvkit_min_target_size ?: 10000
    def cnvkit_method = family_config.cnvkit_method ?: 'cbs'
    def cnvkit_threshold = family_config.cnvkit_threshold
    def cnvkit_ploidy = family_config.cnvkit_ploidy ?: 2
    def cnvkit_drop_outliers = family_config.cnvkit_drop_outliers ?: true
    def expansionhunter_min_anchor = family_config.expansionhunter_min_anchor ?: 8
    def expansionhunter_max_irr = family_config.expansionhunter_max_irr ?: 100
    def tiea_min_support = family_config.tiea_min_support ?: 10
    def tiea_min_softclip = family_config.tiea_min_softclip ?: 36
    def tiea_min_mapq = family_config.tiea_min_mapq ?: 20
    def tiea_cluster_window = family_config.tiea_cluster_window ?: 10
    def tiea_threads = family_config.tiea_threads ?: 4
    def whatshap_chromosomes = family_config.whatshap_chromosomes
    def whatshap_ignore_rg = family_config.whatshap_ignore_rg ?: true
    def whatshap_ref_conf = family_config.whatshap_ref_conf ?: 20
    def vep_use_pick = family_config.vep_use_pick ?: true
    def vep_use_refseq_only = family_config.vep_use_refseq_only ?: false
    def vep_cache_dir = family_config.vep_cache_dir
    def vep_extra_args = family_config.vep_extra_args
    def vep_mei_upstream = family_config.vep_mei_upstream_distance ?: 5000
    def vep_mei_downstream = family_config.vep_mei_downstream_distance ?: 5000
    def vep_db_dir = family_config.vep_db_dir
    def vep_use_gnomad = family_config.vep_use_gnomad ?: true
    def vep_use_alphamissense = family_config.vep_use_alphamissense ?: true
    def vep_use_evoscore = family_config.vep_use_evoscore ?: true
    def vep_use_pangolin = family_config.vep_use_pangolin ?: true
    def vep_flanking_seq_len = family_config.vep_flanking_seq_len ?: 10
    def vep_use_clinvar = family_config.vep_use_clinvar ?: true
    def vep_use_missense_zscore = family_config.vep_use_missense_zscore ?: true
    def stranger_filter_mode = family_config.stranger_filter_mode ?: 'pathogenic'
    def automap_min_roh = family_config.automap_min_roh_length ?: 500000
    def automap_max_roh = family_config.automap_max_roh_length ?: 100000000
    def automap_min_snps = family_config.automap_min_snps_in_roh ?: 20
    def automap_max_gap = family_config.automap_max_gap_length ?: 500000
    def outdir = family_config.outdir ?: './results'

    // =========================================================================
    // 预定义输出通道
    // =========================================================================
    ch_all_alignments = Channel.empty()
    ch_all_alignment_indices = Channel.empty()
    ch_all_gvcfs = Channel.empty()
    ch_all_gvcf_tbis = Channel.empty()
    ch_all_sample_ids = Channel.empty()
    ch_proband_id = Channel.empty()
    ch_father_id = Channel.empty()
    ch_mother_id = Channel.empty()

    // =========================================================================
    // Step 1: 解析家系成员，分离先证者/父/母
    // =========================================================================
    ch_family_reads.tap { tuple ->
        def role = tuple[1]
        def sample_id = tuple[2]

        if (role == 'proband') {
            ch_proband_id = ch_proband_id.mix(Channel.value(sample_id))
        } else if (role == 'father') {
            ch_father_id = ch_father_id.mix(Channel.value(sample_id))
        } else if (role == 'mother') {
            ch_mother_id = ch_mother_id.mix(Channel.value(sample_id))
        }
        ch_all_sample_ids = ch_all_sample_ids.mix(Channel.value(sample_id))
    }

    // =========================================================================
    // Step 2: 所有成员 - FASTQ 质控 (并行)
    // =========================================================================
    FASTP(ch_family_reads.map { family_id, role, sample_id, reads ->
        def meta = [id: sample_id, family_id: family_id, role: role,
                    read_group: "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA"]
        return tuple(meta, reads)
    })
    .publishDir("${outdir}/qc/fastp", mode: 'copy')

    // 提取处理后的 reads 和 sample_id
    ch_processed_reads = FASTP.out.clean_reads
    ch_processed_sample_ids = FASTP.out.clean_reads.map { meta, reads -> meta.id }

    // =========================================================================
    // Step 3: 所有成员 - 序列比对 (并行)
    // =========================================================================
    mem_threshold = 64L * 1024L * 1024L * 1024L
    max_mem_bytes = (params.max_memory as nextflow.util.MemoryUnit).toBytes()
    has_high_resources = max_mem_bytes >= mem_threshold

    if (use_gpu) {
        log.info "比对模式: Parabricks (GPU) - Trio"
        PB_FQ2BAM(
            ch_processed_reads,
            ch_fasta,
            ch_fasta_fai,
            ch_processed_sample_ids
        )
        .publishDir("${outdir}/alignment", mode: 'copy')
        ch_all_alignments = PB_FQ2BAM.out.bam
        ch_all_alignment_indices = PB_FQ2BAM.out.bai
    } else if (has_high_resources) {
        log.info "比对模式: BWA-MEM2 + MarkDuplicates - Trio"
        BWAMEM2(
            ch_processed_reads,
            ch_fasta,
            ch_fasta_fai,
            ch_processed_sample_ids
        )
        .publishDir("${outdir}/alignment", mode: 'copy')
        MARKDUPLICATES(BWAMEM2.out.bam, BWAMEM2.out.bai, ch_fasta)
        .publishDir("${outdir}/alignment", mode: 'copy')
        SAMTOOLS_INDEX(MARKDUPLICATES.out.alignment)
        .publishDir("${outdir}/alignment", mode: 'copy')
        ch_all_alignments = SAMTOOLS_INDEX.out.alignment
        ch_all_alignment_indices = SAMTOOLS_INDEX.out.index
    } else {
        log.info "比对模式: BWA-MEM + MarkDuplicates - Trio"
        BWAMEM(
            ch_processed_reads,
            ch_fasta,
            ch_fasta_fai,
            ch_processed_sample_ids
        )
        .publishDir("${outdir}/alignment", mode: 'copy')
        MARKDUPLICATES(BWAMEM.out.bam, BWAMEM.out.bai, ch_fasta)
        .publishDir("${outdir}/alignment", mode: 'copy')
        SAMTOOLS_INDEX(MARKDUPLICATES.out.alignment)
        .publishDir("${outdir}/alignment", mode: 'copy')
        ch_all_alignments = SAMTOOLS_INDEX.out.alignment
        ch_all_alignment_indices = SAMTOOLS_INDEX.out.index
    }

    // =========================================================================
    // Step 4: 所有成员 - QC 统计 (并行)
    // =========================================================================
    // 覆盖度统计
    ch_cov_alignment = ch_all_alignments
        .combine(ch_all_alignment_indices)
        .combine(ch_processed_sample_ids)
        .map { bam, bai, sid -> tuple(sid, bam) }

    XAMDST(ch_cov_alignment, ch_target_bed)
    .publishDir("${outdir}/qc/coverage", mode: 'copy')

    XAMDST.mt(ch_cov_alignment, ch_mt_bed)
    .publishDir("${outdir}/qc/coverage", mode: 'copy')

    COLLECTQCMETRICS(ch_all_alignments, ch_all_alignment_indices, ch_fasta, ch_target_bed)
    .publishDir("${outdir}/qc/alignment", mode: 'copy')

    SEX_CHECK_SRY(ch_all_alignments, ch_all_alignment_indices, genome_assembly, sex_check_threshold)
    .publishDir("${outdir}/qc/sex_check", mode: 'copy')

    // =========================================================================
    // Step 5: 所有成员 - DeepVariant gVCF (并行)
    // =========================================================================
    if (use_gpu) {
        log.info "变异检测: Parabricks DeepVariant (GPU) - Trio"
        PB_DEEPVARIANT(
            ch_all_alignments,
            ch_all_alignment_indices,
            ch_fasta,
            ch_fasta_fai,
            ch_fasta_dict,
            ch_target_bed,
            deepvariant_model,
            deepvariant_shards
        )
        .publishDir("${outdir}/variants/snv_indel/raw", mode: 'copy')
        ch_all_gvcfs = PB_DEEPVARIANT.out.gvcf
        ch_all_gvcf_tbis = PB_DEEPVARIANT.out.gvcf_tbi
    } else {
        log.info "变异检测: DeepVariant (CPU) - Trio"
        DEEPVARIANT(
            ch_all_alignments,
            ch_all_alignment_indices,
            ch_fasta,
            ch_fasta_fai,
            ch_fasta_dict,
            ch_target_bed,
            deepvariant_model,
            deepvariant_shards
        )
        .publishDir("${outdir}/variants/snv_indel/raw", mode: 'copy')
        ch_all_gvcfs = DEEPVARIANT.out.gvcf
        ch_all_gvcf_tbis = DEEPVARIANT.out.gvcf_tbi
    }

    // =========================================================================
    // Step 6: 家系 gVCF 合并 (GLnexus Trio)
    // =========================================================================
    ch_family_id = ch_family_reads.map { tuple -> tuple[0] }.distinct()

    GLNEXUS_JOINTCALL_TRIO(
        ch_all_gvcfs.collect(),
        ch_all_gvcf_tbis.collect(),
        ch_ped_file,
        glnexus_config,
        ch_family_id.map { fid -> "${fid}_trio" }
    )
    .publishDir("${outdir}/variants/snv_indel/joint", mode: 'copy')

    ch_joint_vcf = GLNEXUS_JOINTCALL_TRIO.out.vcf
    ch_joint_vcf_tbi = GLNEXUS_JOINTCALL_TRIO.out.vcf_tbi

    // =========================================================================
    // Step 7: 亲缘关系验证 (Peddy)
    // =========================================================================
    PEDDY_FAMILY(
        ch_family_id,
        ch_all_gvcfs.collect(),
        ch_all_gvcf_tbis.collect(),
        ch_processed_sample_ids.collect(),
        ch_father_id.collect(),
        ch_mother_id.collect(),
        ch_proband_id.collect(),
        ch_fasta,
        ch_fasta_fai
    )
    .publishDir("${outdir}/qc/peddy", mode: 'copy')

    // =========================================================================
    // Step 8: WhatsHap Trio 定相
    // =========================================================================
    WHATSHAP_PHASE(
        ch_joint_vcf,
        ch_joint_vcf_tbi,
        ch_all_alignments,
        ch_all_alignment_indices,
        ch_fasta,
        ch_fasta_fai,
        ch_processed_sample_ids,
        whatshap_chromosomes,
        whatshap_ignore_rg,
        whatshap_ref_conf,
        "${outdir}/variants/snv_indel/phased"
    )
    .publishDir("${outdir}/variants/snv_indel/phased", mode: 'copy')

    ch_phased_vcf = WHATSHAP_PHASE.out.vcf
    ch_phased_vcf_tbi = WHATSHAP_PHASE.out.vcf_tbi

    // =========================================================================
    // Step 9: VEP 注释
    // =========================================================================
    def db_prefix = genome_assembly == 'GRCh37' ? 'hg19' : 'hg38'
    ch_gnomad_vcf = Channel.empty()
    ch_gnomad_tbi = Channel.empty()
    ch_alphamissense_vcf = Channel.empty()
    ch_alphamissense_tbi = Channel.empty()
    ch_evoscore_vcf = Channel.empty()
    ch_evoscore_tbi = Channel.empty()
    ch_pangolin_vcf = Channel.empty()
    ch_pangolin_tbi = Channel.empty()
    ch_clinvar_vcf = Channel.empty()
    ch_clinvar_tbi = Channel.empty()
    ch_missense_bed = Channel.empty()

    if (vep_db_dir) {
        def db_path = "${vep_db_dir}/${db_prefix}"
        if (vep_use_gnomad) {
            ch_gnomad_vcf = Channel.fromPath("${db_path}/${db_prefix}_gnomad.v4.1.filtered.vcf.gz", checkIfExists: true)
            ch_gnomad_tbi = Channel.fromPath("${db_path}/${db_prefix}_gnomad.v4.1.filtered.vcf.gz.tbi", checkIfExists: true)
        }
        if (vep_use_alphamissense) {
            ch_alphamissense_vcf = Channel.fromPath("${db_path}/${db_prefix}_AlphaMissense.v3.vcf.gz", checkIfExists: true)
            ch_alphamissense_tbi = Channel.fromPath("${db_path}/${db_prefix}_AlphaMissense.v3.vcf.gz.tbi", checkIfExists: true)
        }
        if (vep_use_evoscore) {
            ch_evoscore_vcf = Channel.fromPath("${db_path}/${db_prefix}_EVOScore2.vcf.gz", checkIfExists: true)
            ch_evoscore_tbi = Channel.fromPath("${db_path}/${db_prefix}_EVOScore2.vcf.gz.tbi", checkIfExists: true)
        }
        if (vep_use_pangolin) {
            ch_pangolin_vcf = Channel.fromPath("${db_path}/${db_prefix}_pangolin.vcf.gz", checkIfExists: true)
            ch_pangolin_tbi = Channel.fromPath("${db_path}/${db_prefix}_pangolin.vcf.gz.tbi", checkIfExists: true)
        }
        if (vep_use_clinvar) {
            ch_clinvar_vcf = Channel.fromPath("${vep_db_dir}/clinvar/clinvar.vcf.gz", checkIfExists: true)
            ch_clinvar_tbi = Channel.fromPath("${vep_db_dir}/clinvar/clinvar.vcf.gz.tbi", checkIfExists: true)
        }
        if (vep_use_missense_zscore) {
            ch_missense_bed = Channel.fromPath("${db_path}/missenseByTranscript.${db_prefix}.v4.1.bed", checkIfExists: true)
        }
    }

    VEP_ANNOTATE(
        ch_phased_vcf,
        ch_phased_vcf_tbi,
        ch_fasta,
        ch_fasta_fai,
        genome_assembly,
        vep_use_pick,
        vep_use_refseq_only,
        vep_cache_dir,
        vep_extra_args,
        ch_gnomad_vcf, ch_gnomad_tbi,
        ch_alphamissense_vcf, ch_alphamissense_tbi,
        ch_evoscore_vcf, ch_evoscore_tbi,
        ch_pangolin_vcf, ch_pangolin_tbi,
        vep_flanking_seq_len,
        ch_clinvar_vcf, ch_clinvar_tbi,
        ch_missense_bed
    )
    .publishDir("${outdir}/variants/snv_indel/annotated", mode: 'copy')

    ch_vep_vcf = VEP_ANNOTATE.out.vep_vcf
    ch_vep_vcf_tbi = VEP_ANNOTATE.out.vep_vcf_tbi

    // =========================================================================
    // Step 10: 先证者专属分析
    // =========================================================================
    // 分离先证者数据
    ch_proband_alignment = ch_all_alignments
        .combine(ch_processed_sample_ids)
        .filter { bam, sid -> sid == ch_proband_id.collect()[0] }
        .map { bam, sid -> bam }

    ch_proband_alignment_index = ch_all_alignment_indices
        .combine(ch_processed_sample_ids)
        .filter { bai, sid -> sid == ch_proband_id.collect()[0] }
        .map { bai, sid -> bai }

    ch_proband_sex = SEX_CHECK_SRY.out.json
        .combine(ch_proband_id)
        .filter { json, sid -> json.baseName.contains(sid) }
        .map { json, sid ->
            def content = new groovy.json.JsonSlurper().parseText(json.text)
            tuple(sid, content.inferred_sex == 'Male')
        }

    // CNV
    CNVKIT_BATCH(
        ch_proband_alignment,
        ch_proband_alignment_index,
        ch_target_bed,
        cnvkit_annotate,
        cnvkit_split_size,
        ch_cnvkit_access_bed,
        cnvkit_min_target_size,
        ch_fasta,
        ch_cnvkit_reference,
        cnvkit_method,
        cnvkit_threshold,
        cnvkit_ploidy,
        ch_proband_sex.map { sid, male -> male },
        cnvkit_drop_outliers,
        "${outdir}/proband/cnv"
    )
    .publishDir("${outdir}/proband/cnv", mode: 'copy')

    // STR
    ch_proband_sex_str = ch_proband_sex.map { sid, male -> male ? 'male' : 'female' }

    EXPANSIONHUNTER(
        ch_proband_alignment,
        ch_proband_alignment_index,
        ch_fasta,
        ch_fasta_fai,
        ch_expansionhunter_catalog,
        ch_proband_sex_str,
        expansionhunter_min_anchor,
        expansionhunter_max_irr,
        "${outdir}/proband/str/raw"
    )
    .publishDir("${outdir}/proband/str/raw", mode: 'copy')

    STRANGER_ANNOTATE(
        EXPANSIONHUNTER.out.vcf,
        EXPANSIONHUNTER.out.vcf_tbi,
        ch_expansionhunter_catalog,
        genome_assembly,
        "${outdir}/proband/str/raw"
    )
    .publishDir("${outdir}/proband/str/annotated", mode: 'copy')

    STRANGER_FILTER(
        STRANGER_ANNOTATE.out.vcf,
        STRANGER_ANNOTATE.out.vcf_tbi,
        stranger_filter_mode,
        "${outdir}/proband/str/raw"
    )
    .publishDir("${outdir}/proband/str/filtered", mode: 'copy')

    // MEI
    TIEA_WES(
        ch_proband_alignment,
        ch_proband_alignment_index,
        ch_fasta,
        ch_proband_id,
        tiea_min_support,
        tiea_min_softclip,
        tiea_min_mapq,
        tiea_cluster_window,
        tiea_threads
    )
    .publishDir("${outdir}/proband/mei/raw", mode: 'copy')

    VEP_MEI(
        TIEA_WES.out.vcf,
        TIEA_WES.out.vcf_tbi,
        ch_fasta,
        ch_fasta_fai,
        genome_assembly,
        vep_mei_upstream,
        vep_mei_downstream,
        vep_cache_dir,
        vep_extra_args
    )
    .publishDir("${outdir}/proband/mei/annotated", mode: 'copy')

    // 线粒体
    MUTECT2_MT(
        ch_proband_alignment,
        ch_proband_alignment_index,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        genome_assembly
    )
    .publishDir("${outdir}/proband/mt/raw", mode: 'copy')

    VEP_MT(
        MUTECT2_MT.out.vcf,
        MUTECT2_MT.out.vcf_tbi,
        ch_fasta,
        ch_fasta_fai,
        genome_assembly,
        vep_cache_dir,
        vep_extra_args
    )
    .publishDir("${outdir}/proband/mt/annotated", mode: 'copy')

    // BAF
    BCFTOOLS_BAF_MATRIX(
        ch_proband_alignment,
        ch_proband_alignment_index,
        ch_fasta,
        ch_fasta_fai,
        ch_snp_positions,
        baf_max_depth
    )
    .publishDir("${outdir}/proband/baf", mode: 'copy')

    // =========================================================================
    // Step 11: Trio ROH/UPD 检测
    // =========================================================================
    // 分离各成员 VCF
    ch_member_vcf = ch_vep_vcf.combine(ch_processed_sample_ids).map { vcf, sid -> tuple(sid, vcf) }
    ch_member_vcf_tbi = ch_vep_vcf_tbi.combine(ch_processed_sample_ids).map { tbi, sid -> tuple(sid, tbi) }

    // 先证者 ROH
    AUTOMAP_ROH(
        ch_member_vcf.filter { sid, vcf -> sid == ch_proband_id.collect()[0] }.map { sid, vcf -> vcf },
        ch_member_vcf_tbi.filter { sid, tbi -> sid == ch_proband_id.collect()[0] }.map { sid, tbi -> tbi },
        ch_fasta,
        ch_fasta_fai,
        ch_target_bed,
        genome_assembly,
        automap_min_roh,
        automap_max_roh,
        automap_min_snps,
        automap_max_gap,
        ch_proband_id,
        "${outdir}/homozygosity"
    )
    .publishDir("${outdir}/homozygosity", mode: 'copy')

    // Trio UPD (有父母数据)
    AUTOMAP_UPD(
        ch_member_vcf.filter { sid, vcf -> sid == ch_proband_id.collect()[0] }.map { sid, vcf -> vcf },
        ch_member_vcf_tbi.filter { sid, tbi -> sid == ch_proband_id.collect()[0] }.map { sid, tbi -> tbi },
        ch_fasta,
        ch_fasta_fai,
        genome_assembly,
        ch_proband_id,
        ch_father_id,
        ch_mother_id,
        ch_member_vcf.filter { sid, vcf -> ch_father_id.collect().contains(sid) }.map { sid, vcf -> vcf },
        ch_member_vcf_tbi.filter { sid, tbi -> ch_father_id.collect().contains(sid) }.map { sid, tbi -> tbi },
        ch_member_vcf.filter { sid, vcf -> ch_mother_id.collect().contains(sid) }.map { sid, vcf -> vcf },
        ch_member_vcf_tbi.filter { sid, tbi -> ch_mother_id.collect().contains(sid) }.map { sid, tbi -> tbi },
        "${outdir}/homozygosity"
    )
    .publishDir("${outdir}/homozygosity", mode: 'copy')

    // =========================================================================
    // Emit Results
    // =========================================================================
    emit:
    joint_vcf              = ch_joint_vcf
    joint_vcf_tbi          = ch_joint_vcf_tbi
    de_novo_vcf            = GLNEXUS_JOINTCALL_TRIO.out.de_novo_vcf
    phased_vcf             = ch_phased_vcf
    vep_vcf                = ch_vep_vcf
    peddy_result           = PEDDY_FAMILY.out.family_result
    cnv_results            = CNVKIT_BATCH.out.seg_call_cns
    str_results            = STRANGER_ANNOTATE.out.vcf
    str_filtered           = STRANGER_FILTER.out.filtered_vcf
    mei_results            = TIEA_WES.out.vcf
    mt_vcf                 = MUTECT2_MT.out.vcf
    baf_matrix             = BCFTOOLS_BAF_MATRIX.out.baf_matrix
    roh_results            = AUTOMAP_ROH.out.roh_bed
    upd_results            = AUTOMAP_UPD.out.upd_bed
}