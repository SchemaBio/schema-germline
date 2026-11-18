version 1.2

import "tasks/fastp.wdl" as FASTP
import "tasks/bwamem.wdl" as BWAMEM
import "tasks/sambamba.wdl" as SAMBAMBA

workflow PreBam {
    input {
        String prefix
        File read_1
        File read_2
        String fasta
        Directory ref_dir
        Int threads
    }

    String ref_fasta_name = basename(fasta)

    # BAM文件生产线
    call FASTP.Fastp as Fastp {
        input:
            prefix = prefix,
            read_1 = read_1,
            read_2 = read_2,
            threads = threads
    }
    call BWAMEM.BwaAlign as BwaAlign {
        input:
            prefix = prefix,
            read_1 = Fastp.clean_read_1,
            read_2 = Fastp.clean_read_2,
            ref_dir = ref_dir,
            ref_fasta_name = ref_fasta_name,
            threads = threads
    }
    call SAMBAMBA.SambambaMarkdup as Markdup {
        input:
            prefix = prefix,
            bam = BwaAlign.out_bam,
            bai = BwaAlign.out_bai,
            threads = threads
    }

    output {
        File out_bam = Markdup.markdup_bam
        File out_bai = Markdup.markdup_bai
        File sambamba_stats = Markdup.sambamba_stats
        File qc_report = Fastp.json_report
    }
}
