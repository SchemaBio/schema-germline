version 1.2

# miniwdl run /home/ubuntu/schema-germline/single.wdl -p /home/ubuntu/schema-germline --cfg /home/ubuntu/schema-germline/conf/local.cfg -i /home/ubuntu/schema-germline/inputs/single.json --dir /mnt/data/output
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

# 输出类型定义
struct PipelineSummary {
    String prefix
    String status
    String pipeline
    String version
    File bam
    File bai
    File bed
    File qc_result
    File vcf_raw
    File snp_indel
    File mt
    File cnv_region
    File cnv_gene
    File cnv_raw
    File str
    File mei
    File roh
}

workflow SingleWES {
    input {
        String prefix
        File read_1
        File read_2
        String fasta
        File bed
        Int flank_size
        String assembly
        File cnvkit_reference
        String cnvkit_segmentation_method
        Int sry_sex_cutoff
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

    # BAM文件生产线
    call FASTP.Fastp as Fastp {
        input:
            prefix = prefix,
            read_1 = read_1,
            read_2 = read_2,
            threads = 16
    }
    call BWAMEM.BwaAlign as BwaAlign {
        input:
            prefix = prefix,
            read_1 = Fastp.clean_read_1,
            read_2 = Fastp.clean_read_2,
            ref_dir = ref_dir,
            ref_fasta_name = ref_fasta_name,
            threads = 32
    }
    call SAMBAMBA.SambambaMarkdup as Markdup {
        input:
            prefix = prefix,
            bam = BwaAlign.out_bam,
            bai = BwaAlign.out_bai,
            threads = 32
    }

    # 废弃的Markdup任务，使用Sambamba替代GATK进行去重
    # call GATK.MarkDuplicates as Markdup {
    #     input:
    #         prefix = prefix,
    #         bam = BwaAlign.out_bam,
    #         bai = BwaAlign.out_bai
    # }

    # 质控线路
    call SAMTOOLS.SamtoolsSexCheck as SamtoolsSexCheck {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            assembly = assembly
    }
    call XAMDST.Xamdst as Xamdst {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            bed = FixBed.fixed_bed,
            threads = 8
    }
    call GERMLINE.CreateMitoBed as CreateMitoBed {
        input:
            prefix = prefix
    }
    String mt_xamdst_prefix = "~{prefix}.mt"
    call XAMDST.Xamdst as MtXamdst {
        input:
            prefix = mt_xamdst_prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            bed = CreateMitoBed.mito_bed,
            threads = 8
    }
    call GATK.CollectQCMetrics as CollectQCMetrics {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            bed = FixBed.fixed_bed,
            fasta = ref_fasta_name,
            threads = 8,
            ref_dir = ref_dir
    }
    call GERMLINE.FingerPrint as FingerPrint {
        input:
            prefix = prefix,
            fasta = ref_fasta_name,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            assembly = assembly,
            ref_dir = ref_dir,
            threads = 4
    }
    call GERMLINE.QCReport as QCReport {
        input:
            prefix = prefix,
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

    # SNP InDel 分析
    call DEEPVARIANT.DeepVariant as DeepVariant {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            fasta = ref_fasta_name,
            bed = FixBed.fixed_bed,
            threads = 16,
            flank_size = 50,
            ref_dir = ref_dir
    }
    call WHATSHAP.Whatshap as Whatshap {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            vcf = DeepVariant.vcf,
            fasta = ref_fasta_name,
            threads = 8,
            ref_dir = ref_dir
    }
    call GATK.LeftAlignAndTrimVariants as LeftAlignAndTrimVariants {
        input:
            prefix = prefix,
            vcf = Whatshap.out_vcf,
            vcf_tbi = Whatshap.out_vcf_tbi,
            fasta = ref_fasta_name,
            ref_dir = ref_dir
    }
    call VEP.VEP as VEP {
        input:
            prefix = prefix,
            vcf = LeftAlignAndTrimVariants.left_vcf,
            cache_dir = cache_dir,
            schema_bundle = schema_bundle,
            threads = 16,
            assembly = assembly,
            fasta = ref_fasta_name,
            clinvar_version = clinvar_version,
            ref_dir = ref_dir
    }
    call GERMLINE.SNPInDelReport as SNPInDelReport {
        input:
            prefix = prefix,
            vep_vcf = VEP.out_vcf,
            sry_file = SamtoolsSexCheck.SRY_count,
            sex_cutoff = sry_sex_cutoff
    }

    # 线粒体分析
    call GATK.MitochondrialMutect2 as MitochondrialMutect2 {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
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
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            threads = 8
    }
    call CNVKIT.CNVKitFix as CNVKitFix {
        input:
            prefix = prefix,
            target_coverage = CNVKitCoverage.target_coverage,
            antitarget_coverage = CNVKitCoverage.antitarget_coverage,
            reference = cnvkit_reference,
            segmentation_method = cnvkit_segmentation_method
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
            vcf = LeftAlignAndTrimVariants.left_vcf,
            assembly = assembly,
            threads = 4
    }
    call GERMLINE.ROHReport as ROHReport {
        input:
            prefix = prefix
            automap_report = AutoMap.combined_tsv,
            assembly = assembly
    }

    # 转座子分析
    call TIEAWES.TIEA_WES as TIEA_WES {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai
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
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            fasta = ref_fasta_name,
            sry_file = SamtoolsSexCheck.SRY_count,
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
                pipeline: "WES_Single",
                version: "v0.0.1",
                bam: Markdup.markdup_bam,
                bai: Markdup.markdup_bai,
                bed: FixBed.fixed_bed,
                qc_result: QCReport.qc_result,
                vcf_raw: LeftAlignAndTrimVariants.left_vcf,
                snp_indel: SNPInDelReport.snp_indel_result,
                mt
                cnv_region
                cnv_gene
                cnv_raw: CNVKitFix.cnvkit_cns,
                str: STRReport.str_result,
                mei: MEIReport.mei_result,
                roh: ROHReport.roh_result
            }
        )
    }
}
