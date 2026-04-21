version 1.2

task SamtoolsSexCheck {
    input {
        String prefix
        File bam
        File bai
        String assembly
    }

    command <<<
        if [ "~{assembly}" == "hg38" ]; then
            samtools view ~{bam} Y:2654896-2655723 | wc -l > ~{prefix}.SRY.count.txt
        elif [ "~{assembly}" == "hg19" ]; then
            samtools view ~{bam} Y:2649520-2650357 | wc -l > ~{prefix}.SRY.count.txt
        else
            echo "Unsupported assembly: ~{assembly}" >&2
            exit 1
        fi
    >>>

    output {
        File SRY_count = "~{prefix}.SRY.count.txt"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/mapping:v1.0.0"
    }
}