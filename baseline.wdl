version 1.2

import "tasks/fastp.wdl" as FASTP
import "tasks/germline.wdl" as GERMLINE
import "tasks/cnvkit.wdl" as CNVKIT
import "tasks/bwamem.wdl" as BWAMEM

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
    Int fastp_threads = 8
    Int bwa_threads = 16
    Int cnvkit_threads = 8

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
            fasta = fasta,
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

        call CNVKIT.CNVKitCoverage as CNVKitCoverage {
            input:
                prefix = "~{prefix}_sample~{i}",
                target_bed = TargetBed.target_bed,
                antitarget_bed = CNVKitAntitarget.antitarget_bed,
                bam = BwaAlign.out_bam,
                bai = BwaAlign.out_bai,
                threads = cnvkit_threads
        }
    }

    call CNVKIT.CNVKitReference as CNVKitReference {
        input:
            prefix = prefix,
            fasta = fasta,
            target_coverages = CNVKitCoverage.target_coverage,
            antitarget_coverages = CNVKitCoverage.antitarget_coverage,
            ref_dir = ref_dir
    }
}
