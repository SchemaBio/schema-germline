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
include { CNVKIT_BATCH } from '../../modules/cnvkit/main'
include { EXPANSIONHUNTER } from '../../modules/expansionhunter/main'
include { TIEA_WES } from '../../modules/tiea-wes/main'
include { WHATSHAP_PHASE } from '../../modules/whatshap/main'
include { VEP_ANNOTATE; VEP_MT; VEP_MEI } from '../../modules/vep/main'
include { STRANGER_ANNOTATE; STRANGER_FILTER } from '../../modules/stranger/main'

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
    // CNVkit 参数
    ch_cnvkit_reference      // CNVkit 基线文件 (reference.cnn)
    ch_cnvkit_access_bed     // CNVkit 可访问区域 BED (可选)
    val cnvkit_annotate      // CNVkit 是否添加基因注释 (boolean)
    val cnvkit_split_size    // CNVkit 分割大小 (integer, 默认 5000)
    val cnvkit_min_target_size // CNVkit 最小目标大小 (integer, 默认 10000)
    val cnvkit_method        // CNVkit 分段方法: 'cbs', 'haar', 'flasso' (默认 cbs)
    val cnvkit_threshold     // CNVkit 分段阈值 (可选)
    val cnvkit_ploidy        // CNVkit 倍性 (integer, 默认 2)
    val cnvkit_drop_outliers // CNVkit 是否丢弃离群值 (boolean, 默认 true)
    val cnvkit_output_dir    // CNVkit 输出目录
    // ExpansionHunter 参数
    ch_expansionhunter_catalog // ExpansionHunter STR 位点定义文件 (JSON)
    val expansionhunter_min_anchor // ExpansionHunter 最小锚定长度 (integer, 默认 8)
    val expansionhunter_max_irr // ExpansionHunter 最大不完美匹配距离 (integer, 默认 100)
    val expansionhunter_output_dir // ExpansionHunter 输出目录
    // TIEA-WES 参数
    val tiea_min_support      // TIEA-WES 最小断点支持 reads 数 (integer, 默认 10)
    val tiea_min_softclip      // TIEA-WES 最小 softclip 长度 (integer, 默认 36)
    val tiea_min_mapq          // TIEA-WES 最小 MAPQ (integer, 默认 20)
    val tiea_cluster_window    // TIEA-WES 断点聚类窗口 bp (integer, 默认 10)
    val tiea_threads           // TIEA-WES BWA 线程数 (integer, 默认 4)
    // WhatsHap 参数
    val whatshap_chromosomes    // WhatsHap 定相染色体 (逗号分隔, 默认全部)
    val whatshap_ignore_rg      // WhatsHap 是否忽略 read groups (boolean, 默认 true)
    val whatshap_ref_conf       // WhatsHap 参考置信度阈值 (integer, 默认 20)
    val whatshap_output_dir     // WhatsHap 输出目录
    // VEP 参数
    val vep_use_pick            // VEP 是否选择最显著注释 (boolean, 默认 true)
    val vep_use_refseq_only     // VEP 是否仅使用 RefSeq 转录本 (boolean, 默认 false)
    val vep_cache_dir           // VEP 缓存目录 (可选，默认使用容器内置缓存)
    val vep_extra_args          // VEP 额外参数 (可选)
    val vep_mei_upstream_distance // VEP MEI 注释上游基因距离阈值 bp (默认 5000)
    val vep_mei_downstream_distance // VEP MEI 注释下游基因距离阈值 bp (默认 5000)
    // Stranger 参数
    val stranger_filter_mode    // Stranger 过滤模式: 'pathogenic', 'borderline', 'all_disease' (默认 pathogenic)

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
    ch_cnv = Channel.empty()
    ch_str = Channel.empty()
    ch_mei = Channel.empty()
    ch_phasing = Channel.empty()
    ch_vep_vcf = Channel.empty()
    ch_vep_vcf_tbi = Channel.empty()
    ch_vep_mt_vcf = Channel.empty()
    ch_vep_mt_vcf_tbi = Channel.empty()
    ch_vep_mei_vcf = Channel.empty()
    ch_vep_mei_vcf_tbi = Channel.empty()

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
    // Step 9: 单倍型定相 (WhatsHap)
    // =========================================================================
    // 定相应在 VEP 注释之前，以便注释定相后的 VCF
    WHATSHAP_PHASE(
        ch_vcf,
        ch_vcf_tbi,
        ch_alignment,
        ch_alignment_index,
        ch_fasta,
        ch_fasta_fai,
        ch_sample_id,
        whatshap_chromosomes,
        whatshap_ignore_rg,
        whatshap_ref_conf,
        whatshap_output_dir
    )

    ch_phasing = ch_phasing.mix(WHATSHAP_PHASE.out.vcf)
    ch_phasing = ch_phasing.mix(WHATSHAP_PHASE.out.vcf_tbi)

    // 更新 VCF 通道为定相后的 VCF，用于后续 VEP 注释
    ch_vcf = WHATSHAP_PHASE.out.vcf
    ch_vcf_tbi = WHATSHAP_PHASE.out.vcf_tbi

    // =========================================================================
    // Step 10: CNV 检测 (CNVkit)
    // =========================================================================
    // 从性别检测结果提取样本性别，用于 CNV 分析
    ch_sex_result = SEX_CHECK_SRY.out.json.map { json_file ->
        def content = new groovy.json.JsonSlurper().parseText(json_file.text)
        def sample_id = json_file.baseName.replaceAll(/\\.sex\\.json$/, '')
        tuple(sample_id, content.inferred_sex == 'Male')
    }

    // 组合 alignment、index、性别信息
    ch_alignment_with_sex = ch_alignment
        .combine(ch_alignment_index)
        .combine(ch_sample_id)
        .combine(ch_sex_result)
        .map { alignment, index, sample_id, sex_sample_id, is_male ->
            // 确保 sample_id 匹配
            assert sample_id == sex_sample_id || sex_sample_id.startsWith(sample_id)
            tuple(alignment, index, is_male)
        }

    // 执行 CNV 检测
    CNVKIT_BATCH(
        ch_alignment_with_sex.map { alignment, index, is_male -> alignment },
        ch_alignment_with_sex.map { alignment, index, is_male -> index },
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
        ch_alignment_with_sex.map { alignment, index, is_male -> is_male },
        cnvkit_drop_outliers,
        cnvkit_output_dir
    )
    .publishDir(cnvkit_output_dir ?: 'NO_OUTPUT', mode: 'copy', enabled: cnvkit_output_dir != 'NO_OUTPUT')

    ch_cnv = ch_cnv.mix(CNVKIT_BATCH.out.target_bed)
    ch_cnv = ch_cnv.mix(CNVKIT_BATCH.out.antitarget_bed)
    ch_cnv = ch_cnv.mix(CNVKIT_BATCH.out.target_coverage)
    ch_cnv = ch_cnv.mix(CNVKIT_BATCH.out.antitarget_coverage)
    ch_cnv = ch_cnv.mix(CNVKIT_BATCH.out.cnr)
    ch_cnv = ch_cnv.mix(CNVKIT_BATCH.out.exon_call_cns)
    ch_cnv = ch_cnv.mix(CNVKIT_BATCH.out.seg_cns)
    ch_cnv = ch_cnv.mix(CNVKIT_BATCH.out.seg_call_cns)

    // =========================================================================
    // Step 11: STR 扩展检测 (ExpansionHunter)
    // =========================================================================
    // 从性别检测结果提取样本性别，用于 STR 分析
    ch_sex_for_str = ch_sex_result.map { sample_id, is_male ->
        is_male ? 'male' : 'female'
    }

    // 组合 alignment、index、性别信息
    ch_alignment_for_str = ch_alignment
        .combine(ch_alignment_index)
        .combine(ch_sex_for_str)

    // 执行 STR 扩展检测
    EXPANSIONHUNTER(
        ch_alignment_for_str.map { alignment, index, sex -> alignment },
        ch_alignment_for_str.map { alignment, index, sex -> index },
        ch_fasta,
        ch_fasta_fai,
        ch_expansionhunter_catalog,
        ch_alignment_for_str.map { alignment, index, sex -> sex },
        expansionhunter_min_anchor,
        expansionhunter_max_irr,
        expansionhunter_output_dir
    )
    .publishDir(expansionhunter_output_dir ?: 'NO_OUTPUT', mode: 'copy', enabled: expansionhunter_output_dir != 'NO_OUTPUT')

    ch_str = ch_str.mix(EXPANSIONHUNTER.out.vcf)
    ch_str = ch_str.mix(EXPANSIONHUNTER.out.vcf_tbi)
    ch_str = ch_str.mix(EXPANSIONHUNTER.out.json)
    ch_str = ch_str.mix(EXPANSIONHUNTER.out.reads_json)

    // =========================================================================
    // Step 11b: STR 注释 (Stranger)
    // =========================================================================
    // Stranger 注释 ExpansionHunter 输出的 VCF，添加疾病信息和致病范围
    STRANGER_ANNOTATE(
        EXPANSIONHUNTER.out.vcf,
        EXPANSIONHUNTER.out.vcf_tbi,
        ch_expansionhunter_catalog,
        genome_assembly,
        expansionhunter_output_dir
    )

    ch_str = ch_str.mix(STRANGER_ANNOTATE.out.vcf)
    ch_str = ch_str.mix(STRANGER_ANNOTATE.out.vcf_tbi)
    ch_str = ch_str.mix(STRANGER_ANNOTATE.out.summary)

    // 筛选致病性 STR 扩展
    STRANGER_FILTER(
        STRANGER_ANNOTATE.out.vcf,
        STRANGER_ANNOTATE.out.vcf_tbi,
        stranger_filter_mode,
        expansionhunter_output_dir
    )

    ch_str = ch_str.mix(STRANGER_FILTER.out.filtered_vcf)
    ch_str = ch_str.mix(STRANGER_FILTER.out.filtered_vcf_tbi)
    ch_str = ch_str.mix(STRANGER_FILTER.out.report)

    // =========================================================================
    // Step 12: MEI 检测
    // =========================================================================
    TIEA_WES(
        ch_alignment,
        ch_alignment_index,
        ch_fasta,
        ch_sample_id,
        tiea_min_support,
        tiea_min_softclip,
        tiea_min_mapq,
        tiea_cluster_window,
        tiea_threads
    )

    ch_mei = ch_mei.mix(TIEA_WES.out.vcf)
    ch_mei = ch_mei.mix(TIEA_WES.out.vcf_tbi)

    // =========================================================================
    // Step 14: 变异注释 (VEP)
    // =========================================================================
    // SNV/Indel 注释
    VEP_ANNOTATE(
        ch_vcf,
        ch_vcf_tbi,
        ch_fasta,
        ch_fasta_fai,
        genome_assembly,
        vep_use_pick,
        vep_use_refseq_only,
        vep_cache_dir,
        vep_extra_args
    )
    ch_vep_vcf = ch_vep_vcf.mix(VEP_ANNOTATE.out.vep_vcf)
    ch_vep_vcf_tbi = ch_vep_vcf_tbi.mix(VEP_ANNOTATE.out.vep_vcf_tbi)

    // 线粒体注释
    VEP_MT(
        MUTECT2_MT.out.vcf,
        MUTECT2_MT.out.vcf_tbi,
        ch_fasta,
        ch_fasta_fai,
        genome_assembly,
        vep_cache_dir,
        vep_extra_args
    )
    ch_vep_mt_vcf = ch_vep_mt_vcf.mix(VEP_MT.out.vep_vcf)
    ch_vep_mt_vcf_tbi = ch_vep_mt_vcf_tbi.mix(VEP_MT.out.vep_vcf_tbi)

    // MEI 注释
    VEP_MEI(
        TIEA_WES.out.vcf,
        TIEA_WES.out.vcf_tbi,
        ch_fasta,
        ch_fasta_fai,
        genome_assembly,
        vep_mei_upstream_distance,
        vep_mei_downstream_distance,
        vep_cache_dir,
        vep_extra_args
    )
    ch_vep_mei_vcf = ch_vep_mei_vcf.mix(VEP_MEI.out.vep_vcf)
    ch_vep_mei_vcf_tbi = ch_vep_mei_vcf_tbi.mix(VEP_MEI.out.vep_vcf_tbi)

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
    cnv_results             = ch_cnv                             // CNV 检测结果 (CNVkit)
    str_results             = ch_str                             // STR 扩展检测结果 (ExpansionHunter + Stranger)
    stranger_vcf             = STRANGER_ANNOTATE.out.vcf          // Stranger 注释后的 STR VCF
    stranger_vcf_tbi         = STRANGER_ANNOTATE.out.vcf_tbi      // Stranger VCF 索引
    stranger_summary         = STRANGER_ANNOTATE.out.summary      // Stranger 注释摘要
    stranger_filtered_vcf    = STRANGER_FILTER.out.filtered_vcf   // Stranger 过滤后的致病性 STR VCF
    stranger_filtered_report = STRANGER_FILTER.out.report         // Stranger 致病性 STR 报告
    mei_results             = ch_mei                             // MEI 检测结果 (TIEA-WES)
    phasing_results         = ch_phasing                         // 单倍型定相结果 (WhatsHap)
    phasing_json            = WHATSHAP_PHASE.out.json            // 定相统计 JSON (WhatsHap)
    vep_vcf                 = ch_vep_vcf                         // VEP 注释后的 SNV/Indel VCF
    vep_vcf_tbi             = ch_vep_vcf_tbi                     // VEP SNV/Indel VCF 索引
    vep_summary             = VEP_ANNOTATE.out.summary           // VEP SNV/Indel 注释摘要
    vep_mt_vcf              = ch_vep_mt_vcf                      // VEP 注释后的线粒体 VCF
    vep_mt_vcf_tbi          = ch_vep_mt_vcf_tbi                  // VEP 线粒体 VCF 索引
    vep_mt_summary          = VEP_MT.out.summary                 // VEP 线粒体注释摘要
    vep_mei_vcf             = ch_vep_mei_vcf                     // VEP 注释后的 MEI VCF
    vep_mei_vcf_tbi         = ch_vep_mei_vcf_tbi                 // VEP MEI VCF 索引
    vep_mei_summary         = VEP_MEI.out.summary                // VEP MEI 注释摘要
}