version 1.2

task FixBed {
    input {
        File bed
    }

    command <<<
        python /opt/schema-germline/scripts/process_bed.py -i ~{bed} -o fixed.bed
    >>>

    output {
        File fixed_bed = "fixed.bed"
    }
}

task TargetBed {
    input {
        String prefix
        File bed
        String assembly
    }

    command <<<
        if [ "~{assembly}" == "GRCh38" ]; then
            excludeble_bed="/opt/schema-germline/assets/Gencode.GRCh38.cnvkit.target.bed"
        elif [ "~{assembly}" == "GRCh37" ]; then
            excludeble_bed="/opt/schema-germline/assets/Gencode.GRCh37.cnvkit.target.bed"
        else
            echo "Unsupported assembly: ~{assembly}" >&2
            exit 1
        fi

        bedtools intersect -b ~{bed} -a ${excludeble_bed} -wa -u > ~{prefix}.target.bed
    >>>

    output {
        File target_bed = "~{prefix}.target.bed"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/germline:v0.0.3"
    }

}
