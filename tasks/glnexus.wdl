version 1.2

task GLNexus {
    input {
        String prefix
        Array[File] gvcf
        Array[File] gvcf_tbi
        Int threads
    }

    Int memory_gb = threads * 2

    command <<<
        mkdir -p inputs_dir
        GVCFS=(~{sep=" " gvcf})
        TBIS=(~{sep=" " gvcf_tbi})
        for i in "${!GVCFS[@]}"; do
            ln -s "${GVCFS[$i]}" inputs_dir/
            ln -s "${TBIS[$i]}" inputs_dir/
        done

        # 注意：这里的 inputs_dir/*.g.vcf.gz 假设你的输入文件以后缀 .g.vcf.gz 结尾
        glnexus_cli --config DeepVariantWES \
            --threads ~{threads} \
            inputs_dir/*.g.vcf.gz > ~{prefix}.glnexus.bcf

        echo "Converting BCF to VCF.gz..."
        bcftools view \
            --threads ~{threads} \
            -O z \
            -o ~{prefix}.glnexus.vcf.gz \
            ~{prefix}.glnexus.bcf

        echo "Indexing VCF..."
        bcftools index -t ~{prefix}.glnexus.vcf.gz
    >>>

    output {
        File out_vcf = "~{prefix}.glnexus.vcf.gz"
        File out_tbi = "~{prefix}.glnexus.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/glnexus:v1.4.1"
    }
}