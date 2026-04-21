version 1.2

task Fastp {
    input {
        String prefix
        File read_1
        File read_2
        Int threads
    }

    Int actual_threads = if threads > 16 then 16 else threads
    Int memory_gb = actual_threads * 2

    command <<<
        fastp \
            -i ~{read_1} \
            -I ~{read_2} \
            -o ~{prefix}_R1.clean.fq.gz \
            -O ~{prefix}_R2.clean.fq.gz \
            -w ~{actual_threads} \
            -j ~{prefix}.fastp_stats.json \
            -h ~{prefix}.fastp_stats.html \
            --detect_adapter_for_pe
    >>>

    output {
        File clean_read_1 = "~{prefix}_R1.clean.fq.gz"
        File clean_read_2 = "~{prefix}_R2.clean.fq.gz"
        File json_report = "~{prefix}.fastp_stats.json"
        File html_report = "~{prefix}.fastp_stats.html"
    }

    runtime {
        cpu: actual_threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/mapping:v1.0.0"
    }
}