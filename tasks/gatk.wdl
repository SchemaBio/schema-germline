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
            -R ~{ref_dir}/~{fasta} \
            -L MT \
            --mitochondria-mode \
            -I ~{bam} \
            -O ~{prefix}.mt.vcf.gz

        gatk FilterMutectCalls \
            -V ~{prefix}.mt.vcf.gz \
            -R ~{ref_dir}/~{fasta} \
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
        docker: "schemabio/gatk:4.6.2.0"
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
        set -euo pipefail
        read_vcf() {
            case "$1" in
                *.gz|*.bgz) gzip -cd "$1" ;;
                *) cat "$1" ;;
            esac
        }
        # Split multiallelics expands one input row into one row per ALT.
        expected_records=$(read_vcf "~{vcf}" | awk -F '\t' '!/^#/ && NF>=8 {n+=split($5,alts,",")} END {print n+0}')

        # Splitting/trimming multiallelics can emit a later allele before an
        # earlier indel. Build the index only after sorting all normalized rows.
        gatk LeftAlignAndTrimVariants \
            -R "~{ref_dir}/~{fasta}" \
            -V "~{vcf}" \
            -O normalized.unsorted.vcf.gz \
            --split-multi-allelics \
            --create-output-variant-index false

        mkdir -p sort-tmp
        gatk SortVcf \
            -I normalized.unsorted.vcf.gz \
            -O normalized.sorted.vcf.gz \
            --CREATE_INDEX true \
            --TMP_DIR sort-tmp

        actual_records=$(read_vcf normalized.sorted.vcf.gz | awk '!/^#/ && NF>=8 {n++} END {print n+0}')
        if [ "$actual_records" -ne "$expected_records" ]; then
            printf 'LeftAlignAndTrimVariants record count mismatch: expected=%s actual=%s\n' "$expected_records" "$actual_records" >&2
            exit 1
        fi
        test -s normalized.sorted.vcf.gz.tbi
        # Publish the original output names only after normalization, sorting,
        # indexing and completeness checks all succeed.
        mv normalized.sorted.vcf.gz "~{prefix}.left.vcf.gz"
        mv normalized.sorted.vcf.gz.tbi "~{prefix}.left.vcf.gz.tbi"
    >>>

    output {
        File left_vcf = "~{prefix}.left.vcf.gz"
        File left_vcf_tbi = "~{prefix}.left.vcf.gz.tbi"
    }

    runtime {
        docker: "schemabio/gatk:4.6.2.0"
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

        gatk MarkDuplicates \
            -I ~{bam} \
            -O ~{prefix}.markdup.bam \
            -M ~{prefix}.markdup.metrics.txt \
            --TMP_DIR ./tmp \
            --CREATE_INDEX true
        
        mv ~{prefix}.markdup.bai ~{prefix}.markdup.bam.bai

    >>>

    output {
        File markdup_bam = "~{prefix}.markdup.bam"
        File markdup_bai = "~{prefix}.markdup.bam.bai"
        File metrics = "~{prefix}.markdup.metrics.txt"
    }

    runtime {
        docker: "schemabio/gatk:4.6.2.0"
    }
}

# Xamdst有Insertsize统计，不再用GATK的
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
        # gatk CollectInsertSizeMetrics \
        #     -I ~{bam} \
        #     -O ~{prefix}.insertsize.txt \
        #     -H ~{prefix}.histogram.pdf
        
        # gatk EstimateLibraryComplexity \
        #     -I ~{bam} \
        #     -O ~{prefix}.complexity_metrics.txt

        gatk CollectAlignmentSummaryMetrics \
            -I ~{bam} \
            -R ~{ref_dir}/~{fasta} \
            -O ~{prefix}.metrics.txt

        gatk BedToIntervalList \
            -I ~{bed} \
            -O ~{prefix}.interval_list \
            -SD ~{ref_dir}/~{fasta}

        gatk CollectHsMetrics \
            -BI ~{prefix}.interval_list \
            -TI ~{prefix}.interval_list \
            -I ~{bam} \
            -O ~{prefix}.hs.txt
    >>>

    output {
        File summary = "~{prefix}.metrics.txt"
        File hs_metric = "~{prefix}.hs.txt"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "schemabio/gatk:4.6.2.0"
    }
}

