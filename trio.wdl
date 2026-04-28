version 1.2

# 只适用于子父母三人家系
# miniwdl run /home/ubuntu/schema-germline/trio.wdl -p /home/ubuntu/schema-germline --cfg /home/ubuntu/schema-germline/conf/local.cfg -i /home/ubuntu/schema-germline/inputs/trio.json --dir /mnt/data/output
import "tasks/fastp.wdl" as FASTP
import "tasks/bwamem.wdl" as BWAMEM
import "tasks/gatk.wdl" as GATK
import "tasks/whatshap.wdl" as WHATSHAP
import "tasks/vep.wdl" as VEP
import "tasks/sambamba.wdl" as SAMBAMBA
import "tasks/samtools.wdl" as SAMTOOLS
import "tasks/xamdst.wdl" as XAMDST
import "tasks/deepvariant.wdl" as DEEPVARIANT
import "tasks/germline.wdl" as GERMLINE
import "tasks/tiea_wes.wdl" as TIEAWES
import "tasks/cnvkit.wdl" as CNVKIT
import "tasks/expansionhunter.wdl" as EXPANSIONHUNTER
import "tasks/stranger.wdl" as STRANGER
import "tasks/automap.wdl" as AUTOMAP
import "tasks/cnvanno.wdl" as CNVANNO
import "tasks/glnexus.wdl" as GLNEXUS
import "tasks/peddy.wdl" as PEDDY

struct MetaFasta {
    Array[String] members
    Array[File] read_1
    Array[File] read_2
}

# 输出类型定义
struct PipelineSummary {
    String prefix
    String status
    String pipeline
    String version
    Array[File] bam
    Array[File] bai
    File bed
    Array[File] qc_result
    File vcf_raw
    File snp_indel
    File mt
    File cnv_region
    File cnv_gene
    File cnv_raw
    File str
    File mei
    File roh
    File upd
    File peddy
}

workflow TrioWES {
    input {
        String prefix
        MetaFasta meta_info
        String fasta
        File bed
        Int flank_size
        String assembly
        File cnvkit_reference
        File ped
        Int sry_sex_cutoff
        Float cnv_del_threshold = 1.5
        Float cnv_dup_threshold = 2.5
        Int cnv_bin_size = 20000
        Directory ref_dir
        Directory cache_dir
        Directory schema_bundle
    }

    # 参数调整
    String ref_fasta_name = basename(fasta)
    String clinvar_version = '20260415'

    # bed 文件调整
    call GERMLINE.FixBed as FixBed {
        input:
            bed = bed
    }
    call GERMLINE.TargetBed as TargetBed {
        input:
            prefix = prefix,
            bed = FixBed.fixed_bed,
            assembly = assembly
    }
    call GERMLINE.CreateMitoBed as CreateMitoBed {
        input:
            prefix = prefix
    }

    # 批量管线
    scatter (i in range(length(meta_info.members))) {

        # BAM文件生产线
        call FASTP.Fastp as Fastp {
            input:
                prefix = meta_info.members[i],
                read_1 = meta_info.read_1[i],
                read_2 = meta_info.read_2[i],
                threads = 16
        }
        call BWAMEM.BwaAlign as BwaAlign {
            input:
                prefix = meta_info.members[i],
                read_1 = Fastp.clean_read_1,
                read_2 = Fastp.clean_read_2,
                ref_dir = ref_dir,
                ref_fasta_name = ref_fasta_name,
                threads = 16
        }
        call SAMBAMBA.SambambaMarkdup as Markdup {
            input:
                prefix = meta_info.members[i],
                bam = BwaAlign.out_bam,
                bai = BwaAlign.out_bai,
                threads = 16
        }

        # 质控线路
        call SAMTOOLS.SamtoolsSexCheck as SamtoolsSexCheck {
            input:
                prefix = meta_info.members[i],
                bam = Markdup.markdup_bam,
                bai = Markdup.markdup_bai,
                assembly = assembly
        }
        call XAMDST.Xamdst as Xamdst {
            input:
                prefix = meta_info.members[i],
                bam = Markdup.markdup_bam,
                bai = Markdup.markdup_bai,
                bed = FixBed.fixed_bed,
                threads = 8
        }
        call XAMDST.Xamdst as MtXamdst {
            input:
                prefix = meta_info.members[i],
                bam = Markdup.markdup_bam,
                bai = Markdup.markdup_bai,
                bed = CreateMitoBed.mito_bed,
                threads = 8
        }
        call GATK.CollectQCMetrics as CollectQCMetrics {
            input:
                prefix = meta_info.members[i],
                bam = Markdup.markdup_bam,
                bai = Markdup.markdup_bai,
                bed = FixBed.fixed_bed,
                fasta = ref_fasta_name,
                threads = 8,
                ref_dir = ref_dir
        }
        call GERMLINE.FingerPrint as FingerPrint {
            input:
                prefix = meta_info.members[i],
                fasta = ref_fasta_name,
                bam = Markdup.markdup_bam,
                bai = Markdup.markdup_bai,
                assembly = assembly,
                ref_dir = ref_dir,
                threads = 4
        }
        call GERMLINE.QCReport as QCReport {
            input:
                prefix = meta_info.members[i],
                fastp_stats = Fastp.json_report,
                xamdst_json = Xamdst.xamdst_report,
                mt_xamdst_json = MtXamdst.xamdst_report,
                fingerprint_result = FingerPrint.fingerprint_json,
                gatk_metric = CollectQCMetrics.summary,
                hs_metric = CollectQCMetrics.hs_metric,
                sry_result = SamtoolsSexCheck.SRY_count,
                sry_cutoff = sry_sex_cutoff,
                sambamba_stats = Markdup.sambamba_stats
        }

        # DeepVariant 变异检测
        call DEEPVARIANT.DeepVariant as DeepVariant {
            input:
                prefix = meta_info.members[i],
                bam = Markdup.markdup_bam,
                bai = Markdup.markdup_bai,
                fasta = ref_fasta_name,
                bed = FixBed.fixed_bed,
                threads = 16,
                flank_size = flank_size,
                ref_dir = ref_dir
        }
    }

    # SNP InDel 分析 - GLNexus 合并各样本 GVCF
    call GLNEXUS.GLNexus as GLNexus {
        input:
            prefix = prefix,
            gvcf = DeepVariant.gvcf,
            gvcf_tbi = DeepVariant.gvcf_tbi,
            threads = 16
    }

    # 亲缘关系计算
    call PEDDY.Peddy as Peddy {
        input:
            prefix = prefix,
            vcf = GLNexus.out_vcf,
            vcf_idx = GLNexus.out_tbi,
            ped = ped,
            assembly = assembly,
            threads = 8
    }

    # 多线程的whatshap加速处理
    call GERMLINE.SplitVcf as SplitVcfHap {
        input:
            vcf = GLNexus.out_vcf
    }

    scatter (i in range(length(SplitVcfHap.split_vcfs))) {
        call WHATSHAP.Whatshap as Whatshap {
            input:
                prefix = "~{prefix}.part~{i}",
                bam = Markdup.markdup_bam[0],
                bai = Markdup.markdup_bai[0],
                vcf = SplitVcfHap.split_vcfs[i],
                vcf_tbi = SplitVcfHap.split_vcf_tbis[i],
                fasta = ref_fasta_name,
                threads = 2,
                ref_dir = ref_dir
        }
    }

    String phase_vcf_prefix = "~{prefix}.phase.merged"
    call GERMLINE.UniversalMergeVcfs as UniversalMergeVcfsHap {
        input:
            prefix = phase_vcf_prefix,
            vcfs = Whatshap.out_vcf,
            threads = 8
    }

    call GATK.LeftAlignAndTrimVariants as LeftAlignAndTrimVariants {
        input:
            prefix = prefix,
            vcf = UniversalMergeVcfsHap.merged_vcf,
            vcf_tbi = UniversalMergeVcfsHap.merged_vcf_tbi,
            fasta = ref_fasta_name,
            ref_dir = ref_dir
    }

    # 多线程的VEP加速处理
    call GERMLINE.SplitVcf as SplitVcf {
        input:
            vcf = LeftAlignAndTrimVariants.left_vcf
    }

    scatter (i in range(length(SplitVcf.split_vcfs))) {
        call VEP.VEP as VEP_Parallel {
            input:
                prefix = "~{prefix}.part~{i}",
                vcf = SplitVcf.split_vcfs[i],
                cache_dir = cache_dir,
                schema_bundle = schema_bundle,
                threads = 2,
                assembly = assembly,
                fasta = ref_fasta_name,
                clinvar_version = clinvar_version,
                ref_dir = ref_dir
        }
    }

    String merged_vcf_prefix = "~{prefix}.vep.merged"
    call GERMLINE.UniversalMergeVcfs as UniversalMergeVcfs {
        input:
            prefix = merged_vcf_prefix,
            vcfs = VEP_Parallel.out_vcf,
            threads = 8
    }

    call GERMLINE.SNPInDelReport as SNPInDelReport {
        input:
            prefix = prefix,
            vep_vcf = UniversalMergeVcfs.merged_vcf,
            sry_file = SamtoolsSexCheck.SRY_count[0],
            sex_cutoff = sry_sex_cutoff
    }

    # 线粒体分析
    call GATK.MitochondrialMutect2 as MitochondrialMutect2 {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam[0],
            bai = Markdup.markdup_bai[0],
            fasta = ref_fasta_name,
            ref_dir = ref_dir
    }
    String mt_vep_prefix = "~{prefix}.mt"
    call VEP.VEP as MtVEP {
        input:
            prefix = mt_vep_prefix,
            vcf = MitochondrialMutect2.vcf,
            cache_dir = cache_dir,
            schema_bundle = schema_bundle,
            threads = 8,
            assembly = assembly,
            fasta = ref_fasta_name,
            clinvar_version = clinvar_version,
            ref_dir = ref_dir
    }
    call GERMLINE.MTReport as MTReport {
        input:
            prefix = prefix,
            mt_vep_vcf = MtVEP.out_vcf
    }

    # CNV分析
    call CNVKIT.CNVKitAntitarget as CNVKitAntitarget {
        input:
            prefix = prefix,
            fasta = ref_fasta_name,
            assembly = assembly,
            target_bed = TargetBed.target_bed,
            ref_dir = ref_dir
    }
    call CNVKIT.CNVKitCoverage as CNVKitCoverage {
        input:
            prefix = prefix,
            target_bed = TargetBed.target_bed,
            antitarget_bed = CNVKitAntitarget.antitarget_bed,
            bam = Markdup.markdup_bam[0],
            bai = Markdup.markdup_bai[0],
            threads = 8
    }
    call CNVKIT.CNVKitFix as CNVKitFix {
        input:
            prefix = prefix,
            target_coverage = CNVKitCoverage.target_coverage,
            antitarget_coverage = CNVKitCoverage.antitarget_coverage,
            reference = cnvkit_reference
    }
    call GERMLINE.CNVGene as CNVGene {
        input:
            prefix = prefix,
            cnv_cnr = CNVKitFix.cnvkit_cnr,
            cnv_del_threshold = cnv_del_threshold,
            cnv_dup_threshold = cnv_dup_threshold
    }
    call GERMLINE.CNVRegion as CNVRegion {
        input:
            prefix = prefix,
            cnv_cnr = CNVKitFix.cnvkit_cnr,
            cnv_del_threshold = cnv_del_threshold,
            cnv_dup_threshold = cnv_dup_threshold,
            bin_size = cnv_bin_size
    }
    String cnv_gene_prefix = "~{prefix}.gene"
    call CNVANNO.CNVAnno as CNVAnnoGene {
        input:
            prefix = cnv_gene_prefix,
            cnv_bed = CNVGene.cnv_result,
            assembly = assembly
    }
    String cnv_region_prefix = "~{prefix}.region"
    call CNVANNO.CNVAnno as CNVAnnoRegion {
        input:
            prefix = cnv_region_prefix,
            cnv_bed = CNVRegion.cnv_result,
            assembly = assembly
    }

    # LOH：杂合性缺失
    # ROH：连续纯合区域，可能提示隐性遗传病位点
    # UPD：单亲二倍体，可能提示隐性遗传病位点或印记基因异常
    # UPD是ROH的极端情况，ROH是UPD的轻度情况
    # ROH可以由父母双方遗传的相同片段组成，也可以由单亲二倍体（UPD）组成。UPD会导致整个染色体或大片段的ROH
    # 而普通ROH通常较短且分布在多个染色体上。ROH和UPD都可能提示隐性遗传病位点
    # 但UPD还可能涉及印记基因异常。ROH分析可以帮助识别潜在的隐性遗传病位点，而UPD分析可以进一步揭示是否存在单亲二倍体现象
    
    call AUTOMAP.AutoMap as AutoMap {
        input:
            prefix = prefix,
            vcf = DeepVariant.vcf[0],
            assembly = assembly,
            threads = 4
    }
    call GERMLINE.ROHReport as ROHReport {
        input:
            prefix = prefix,
            automap_report = AutoMap.combined_tsv,
            assembly = assembly
    }
    call AUTOMAP.AutoMap as AutoMapParent1 {
        input:
            prefix = meta_info.members[1],
            vcf = DeepVariant.vcf[1],
            assembly = assembly,
            threads = 4
    }
    call GERMLINE.ROHReport as ROHReportParent1 {
        input:
            prefix = meta_info.members[1],
            automap_report = AutoMapParent1.combined_tsv,
            assembly = assembly
    }
    call AUTOMAP.AutoMap as AutoMapParent2 {
        input:
            prefix = meta_info.members[2],
            vcf = DeepVariant.vcf[2],
            assembly = assembly,
            threads = 4
    }
    call GERMLINE.ROHReport as ROHReportParent2 {
        input:
            prefix = meta_info.members[2],
            automap_report = AutoMapParent2.combined_tsv,
            assembly = assembly
    }

    # UPD
    call GERMLINE.UPD as UPD {
        input:
            prefix = prefix,
            proband_roh_report = ROHReport.roh_result,
            father_roh_report = ROHReportParent1.roh_result,
            mother_roh_report = ROHReportParent2.roh_result
    }


    # 转座子分析
    call TIEAWES.TIEA_WES as TIEA_WES {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam[0],
            bai = Markdup.markdup_bai[0]
    }
    String mei_vep_prefix = "~{prefix}.mei"
    call VEP.VEP as MeiVEP {
        input:
            prefix = mei_vep_prefix,
            vcf = TIEA_WES.out_vcf,
            cache_dir = cache_dir,
            schema_bundle = schema_bundle,
            threads = 8,
            assembly = assembly,
            fasta = ref_fasta_name,
            clinvar_version = clinvar_version,
            ref_dir = ref_dir
    }
    call GERMLINE.MEIReport as MEIReport {
        input:
            prefix = prefix,
            mei_vcf = MeiVEP.out_vcf
    }

    # STR分析
    String str_prefix = "~{prefix}.str"
    call EXPANSIONHUNTER.ExpansionHunter as ExpansionHunter {
        input:
            prefix = str_prefix,
            bam = Markdup.markdup_bam[0],
            bai = Markdup.markdup_bai[0],
            fasta = ref_fasta_name,
            sry_file = SamtoolsSexCheck.SRY_count[0],
            ref_dir = ref_dir,
            assembly = assembly,
            threads = 8,
            sex_cutoff = sry_sex_cutoff
    }
    call STRANGER.Stranger as Stranger {
        input:
            prefix = prefix,
            vcf = ExpansionHunter.str_vcf,
            assembly = assembly
    }
    call GERMLINE.STRReport as STRReport {
        input:
            prefix = prefix,
            str_vcf = Stranger.anno_vcf
    }

    output {
        File summary = write_json(
            PipelineSummary {
                prefix: prefix,
                status: "Success",
                pipeline: "WES_Trio",
                version: "v0.0.1",
                bam: Markdup.markdup_bam,
                bai: Markdup.markdup_bai,
                bed: FixBed.fixed_bed,
                qc_result: QCReport.qc_result,
                vcf_raw: LeftAlignAndTrimVariants.left_vcf,
                snp_indel: SNPInDelReport.snp_indel_result,
                mt: MTReport.mt_result,
                cnv_region: CNVAnnoRegion.cnv_anno_result,
                cnv_gene: CNVAnnoGene.cnv_anno_result,
                cnv_raw: CNVKitFix.cnvkit_cnr,
                str: STRReport.str_result,
                mei: MEIReport.mei_result,
                roh: ROHReport.roh_result,
                upd: UPD.upd_result,
                peddy: Peddy.ped_check
            }
        )
    }
}
