version 1.1

task Xamdst {
    input {
        String prefix
        File bam
        File bai
        File bed
        Int threads
    }

    Int memory_gb = threads * 2

    command <<<
        xamdst -1 -p ~{bed} -o result --threads ~{threads} --cutoffdepth 1000 ~{bam}
        cp result/coverage.report ~{prefix}.xamdst.report.json
    >>>

    output {
        File xamdst_report = "~{prefix}.xamdst.report.json"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/mapping:v1.0.0"
    }
}