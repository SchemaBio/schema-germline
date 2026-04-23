version 1.2

task TIEA_WES {
    input {
        String prefix
        File bam
        File bai
    }

    command <<<
        python /app/TIEA-WES.py -p ~{prefix} -i ~{bam} -o result
        cp result/~{prefix}.te.result.vcf ~{prefix}.mei.vcf
    >>>

    output {
        File out_vcf = "~{prefix}.mei.vcf"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/tiea_wes:2.0.1"
    }
}
