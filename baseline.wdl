version 1.2

# miniwdl run /home/ubuntu/schema-germline/baseline.wdl -p /home/ubuntu/schema-germline --cfg /home/ubuntu/schema-germline/conf/local.cfg -i /home/ubuntu/schema-germline/inputs/baseline.json --dir /mnt/data/output
import "tasks/fastp.wdl" as FASTP
import "tasks/germline.wdl" as GERMLINE
import "tasks/cnvkit.wdl" as CNVKIT
import "tasks/bwamem.wdl" as BWAMEM
import "tasks/sambamba.wdl" as SAMBAMBA

workflow CNVBaseline {
    input {
        String prefix
        File bed
        String fasta
        String assembly
        Array[File] read_1
        Array[File] read_2
        Directory ref_dir
    }

    String ref_fasta_name = basename(fasta)
    Int total_threads = 32
    Int num_samples = length(read_1)

    # 动态计算基础线程数
    Int base_threads = total_threads / num_samples

    # 应用各任务的最小/最大约束
    Int fastp_threads = if base_threads < 8 then 8 else if base_threads > 16 then 16 else base_threads
    Int bwa_threads = if base_threads < 16 then 16 else if base_threads > 32 then 32 else base_threads
    Int cnvkit_threads = if base_threads < 8 then 8 else if base_threads > 16 then 16 else base_threads

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

    call CNVKIT.CNVKitAntitarget as CNVKitAntitarget {
        input:
            prefix = prefix,
            target_bed = TargetBed.target_bed,
            fasta = ref_fasta_name,
            assembly = assembly,
            ref_dir = ref_dir
    }

    scatter (i in range(length(read_1))) {
        call FASTP.Fastp as Fastp {
            input:
                prefix = "~{prefix}_sample~{i}",
                read_1 = read_1[i],
                read_2 = read_2[i],
                threads = fastp_threads
        }

        call BWAMEM.BwaAlign as BwaAlign {
            input:
                prefix = "~{prefix}_sample~{i}",
                read_1 = Fastp.clean_read_1,
                read_2 = Fastp.clean_read_2,
                ref_dir = ref_dir,
                ref_fasta_name = ref_fasta_name,
                threads = bwa_threads
        }

        call SAMBAMBA.SambambaMarkdup as Markdup {
            input:
                prefix = "~{prefix}_sample~{i}",
                bam = BwaAlign.out_bam,
                bai = BwaAlign.out_bai,
                threads = bwa_threads
        }
        
        call CNVKIT.CNVKitCoverage as CNVKitCoverage {
            input:
                prefix = "~{prefix}_sample~{i}",
                target_bed = TargetBed.target_bed,
                antitarget_bed = CNVKitAntitarget.antitarget_bed,
                bam = Markdup.markdup_bam,
                bai = Markdup.markdup_bai,
                threads = cnvkit_threads
        }
    }

    call CNVKIT.CNVKitReference as CNVKitReference {
        input:
            prefix = prefix,
            fasta = ref_fasta_name,
            target_coverages = CNVKitCoverage.target_coverage,
            antitarget_coverages = CNVKitCoverage.antitarget_coverage,
            ref_dir = ref_dir
    }
}
