version 1.2

task DeepVariant {
    input {
        String prefix
        File bam
        File bai
        String fasta
        Int threads
        File bed
        Int flank_size = 50
        Directory ref_dir
    }

    Int memory_gb = threads * 2

    command <<<
        awk 'BEGIN {OFS="\t"} {start=$2; end=$3; new_start=start-~{flank_size}; new_end=end+~{flank_size}; if(new_start<0) new_start=0; print $1,new_start,new_end}' ~{bed} > extended.bed

        /opt/deepvariant/bin/run_deepvariant \
            --model_type WES \
            --ref ~{fasta} \
            --reads ~{bam} \
            --output_vcf ~{prefix}.vcf.gz \
            --output_gvcf ~{prefix}.gvcf.vcf.gz \
            --num_shards ~{threads} \
            --regions extended.bed
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File vcf_tbi = "~{prefix}.vcf.gz.tbi"
        File gvcf = "~{prefix}.g.vcf.gz"
        File gvcf_tbi = "~{prefix}.g.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/deepvariant:1.10.0"
    }
}