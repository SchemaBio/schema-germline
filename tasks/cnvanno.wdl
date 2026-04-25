version 1.2

task CNVAnno {
    input {
        String prefix
        File cnv_bed
        String assembly
    }

    command <<<        
        python /app/CNVAnno/cnvanno.py \
            ~{cnv_bed} \
            -d /app/CNVAnno/data \
            -g ~{assembly} \
            -f tsv \
            -o ~{prefix}.cnvanno.txt
    >>>

    output {
        File cnv_anno_result = "~{prefix}.cnvanno.txt"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/cnvanno:v0.0.2"
    }
}