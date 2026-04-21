version 1.2

task MitochondrialMutect2 {
    input {
        String prefix
        File bam
        File bai
        String fasta
        Directory ref_dir
    }

    command <<<
        gatk Mutect2 \
            -R ~{fasta} \
            -L MT \
            --mitochondria-mode \
            -I ~{bam} \
            -O ~{prefix}.mt.vcf.gz

        gatk FilterMutectCalls \
            -V ~{prefix}.mt.vcf.gz \
            -R ~{fasta} \
            -O ~{prefix}.mt.filtered.vcf.gz
        
        gatk SelectVariants \
            -V ~{prefix}.mt.filtered.vcf.gz \
            --exclude-filtered true \
            -O ~{prefix}.mt.pass.vcf.gz
    >>>

    output {
        File vcf = "~{prefix}.mt.vcf.gz"
        File vcf_tbi = "~{prefix}.mt.vcf.gz.tbi"
        File pass_vcf = "~{prefix}.mt.pass.vcf.gz"
        File pass_vcf_tbi = "~{prefix}.mt.pass.vcf.gz.tbi"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/gatk:4.6.2.0"
    }
}

task LeftAlignAndTrimVariants {
    input {
        String prefix
        File vcf
        File vcf_tbi
        String fasta
        Directory ref_dir
    }

    command <<<
        gatk LeftAlignAndTrimVariants \
            -R ~{fasta} \
            -V ~{vcf} \
            -O ~{prefix}.left.vcf.gz \
            --split-multi-allelics
    >>>

    output {
        File left_vcf = "~{prefix}.left.vcf.gz"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/gatk:4.6.2.0"
    }
}

task MarkDuplicates {
    input {
        String prefix
        File bam
        File bai
    }

    command <<<
        set -e
        export TMPDIR=~{tmp_dir}

        gatk MarkDuplicates \
            -I ~{bam} \
            -O ~{prefix}.markdup.bam \
            -M ~{prefix}.markdup.metrics.txt \
            --TMP_DIR ~{tmp_dir} \
            --CREATE_INDEX true
        
        mv ~{prefix}.markdup.bai ~{prefix}.markdup.bam.bai

    >>>

    output {
        File markdup_bam = "~{prefix}.markdup.bam"
        File markdup_bai = "~{prefix}.markdup.bam.bai"
        File metrics = "~{prefix}.markdup.metrics.txt"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/gatk:4.6.2.0"
    }
}

task CollectQCMetrics {
    input {
        String prefix
        File bam
        File bai
        File bed
        String fasta
        Int threads
        Directory ref_dir
    }

    Int memory_gb = threads * 2

    command <<<
        gatk CollectInsertSizeMetrics \
            -I ~{bam} \
            -O ~{prefix}.insertsize.txt \
            -H ~{prefix}.histogram.pdf

        gatk CollectAlignmentSummaryMetrics \
            -I ~{bam} \
            -R ~{fasta} \
            -O ~{prefix}.metrics.txt

        gatk BedToIntervalList \
            -I ~{bed} \
            -O ~{prefix}.interval_list \
            -SD ~{fasta}

        gatk CollectHsMetrics \
            -BI ~{prefix}.interval_list \
            -TI ~{prefix}.interval_list \
            -I ~{bam} \
            -O ~{prefix}.hs.txt

    >>>

    output {
        File insertsizeFile = "~{prefix}.insertsize.txt"
        File summary = "~{prefix}.metrics.txt"
        File hs_metric = "~{prefix}.hs.txt"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/gatk:4.6.2.0"
    }
}

