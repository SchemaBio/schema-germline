version 1.2

task Whatshap {
    input {
        String prefix
        File bam
        File bai
        File vcf
        File vcf_tbi
        String fasta
        Int threads
        Directory ref_dir
    }

    Int memory_gb = threads * 2

    command <<<
        whatshap phase \
            --reference=~{ref_dir}/~{fasta} \
            -o ~{prefix}.phase.vcf \
            ~{vcf} \
            ~{bam}
        bgzip -f ~{prefix}.phase.vcf
        tabix -f -p vcf ~{prefix}.phase.vcf.gz
    >>>

    output {
        File out_vcf = "~{prefix}.phase.vcf.gz"
        File out_vcf_tbi = "~{prefix}.phase.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/whatshap:2.8"
    }
}