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

workflow SingleWES {
    input {
        String prefix
        File read_1
        File read_2
        String fasta
        File bed
        Int flank_size
        String assembly
        Directory ref_dir
        Directory cache_dir
        Directory schema_bundle
        Directory tmp_dir
    }

    # 参数调整
    String ref_fasta_name = basename(fasta)
    File mito_bed = "assets/mito.bed"
    String clinvar_version = '20260415'

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
            bed = bed,
            threads = 8
    }
    call XAMDST.Xamdst as MtXamdst {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            bed = mito_bed,
            threads = 8
    }
    call GATK.CollectQCMetrics as CollectQCMetrics {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            bed = bed,
            fasta = fasta,
            threads = 8,
            ref_dir = ref_dir
    }

    # SNP指纹


    # SNP InDel 分析
    call DEEPVARIANT.DeepVariant as DeepVariant {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            fasta = fasta,
            bed = bed,
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
            fasta = fasta,
            threads = 8,
            ref_dir = ref_dir
    }
    call GATK.LeftAlignAndTrimVariants as LeftAlignAndTrimVariants {
        input:
            prefix = prefix,
            vcf = Whatshap.out_vcf,
            vcf_tbi = Whatshap.out_vcf_tbi,
            fasta = fasta,
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
            fasta = fasta,
            clinvar_version = clinvar_version,
            ref_dir = ref_dir
    }

    # 线粒体分析
    call GATK.MitochondrialMutect2 as MitochondrialMutect2 {
        input:
            prefix = prefix,
            bam = Markdup.markdup_bam,
            bai = Markdup.markdup_bai,
            fasta = fasta,
            ref_dir = ref_dir
    }
    call VEP.VEP as MtVEP {
        input:
            prefix = prefix,
            vcf = MitochondrialMutect2.vcf,
            cache_dir = cache_dir,
            schema_bundle = schema_bundle,
            threads = 16,
            assembly = assembly,
            fasta = fasta,
            clinvar_version = clinvar_version,
            ref_dir = ref_dir
    }

    # STR分析




    # CNV分析



    # ROH / UPD 分析




    # 转座子分析






}
