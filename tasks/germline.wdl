version 1.2

task FixBed {
    input {
        File bed
    }

    command <<<
        python /pipeline/schema-germline/scripts/process_bed.py -i ~{bed} -o fixed.bed
    >>>

    output {
        File fixed_bed = "fixed.bed"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/germline:v1.0.0"
    }
}
