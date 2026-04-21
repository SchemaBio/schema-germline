version 1.2

task SambambaMarkdup {
    input {
        String prefix
        File bam
        File bai
        Int threads
    }

    Int memory_gb = threads * 2

    command <<<
        set -e
        export TMPDIR=~{tmp_dir}

        sambamba markdup \
            --nthreads ~{threads} \
            --overflow-list-size 1000000 \
            --hash-table-size 1000000 \
            --compression-level 1 \
            ~{bam} \
            ~{prefix}.markdup.bam

        samtools index ~{prefix}.markdup.bam
    >>>

    output {
        File markdup_bam = "~{prefix}.markdup.bam"
        File markdup_bai = "~{prefix}.markdup.bam.bai"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/mapping:v1.0.0"
    }
}

