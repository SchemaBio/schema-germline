version 1.2

# DeepVariant的家系版本，显著提升了家系数据的变异检测能力，尤其是de novo变异的检测能力。适用于父母子样本
# 输入的Array[File]中严格落实子父母排序
task DeepTrio {
    input {
        MetaFasta meta_info
        Array[File] bam
        Array[File] bai
        String fasta
        Int threads
        File bed
        Int flank_size = 50
        Directory ref_dir
    }

    Int memory_gb = threads * 2

    command <<<
        awk 'BEGIN {OFS="\t"} {start=$2; end=$3; new_start=start-~{flank_size}; new_end=end+~{flank_size}; if(new_start<0) new_start=0; print $1,new_start,new_end}' ~{bed} > extended.bed

        /opt/deepvariant/bin/deeptrio/run_deeptrio \
            --model_type=WES \
            --ref=~{ref_dir}/~{fasta} \
            --reads_child=~{bam[0]} \
            --reads_parent1=~{bam[1]} \
            --reads_parent2=~{bam[2]} \
            --output_vcf_child ~{meta_info.members[0]}.vcf.gz \
            --output_vcf_parent1 ~{meta_info.members[1]}.vcf.gz \
            --output_vcf_parent2 ~{meta_info.members[2]}.vcf.gz \
            --sample_name_child '~{meta_info.members[0]}' \
            --sample_name_parent1 '~{meta_info.members[1]}' \
            --sample_name_parent2 '~{meta_info.members[2]}' \
            --num_shards ~{threads}  \
            --regions extended.bed \
            --intermediate_results_dir ./intermediate \
            --output_gvcf_child ~{meta_info.members[0]}.g.vcf.gz \
            --output_gvcf_parent1 ~{meta_info.members[1]}.g.vcf.gz \
            --output_gvcf_parent2 ~{meta_info.members[2]}.g.vcf.gz
    >>>

    output {
        File child_vcf = "~{meta_info.members[0]}.vcf.gz"
        File child_vcf_index = "~{meta_info.members[0]}.vcf.gz.tbi"
        File parent1_vcf = "~{meta_info.members[1]}.vcf.gz"
        File parent1_vcf_index = "~{meta_info.members[1]}.vcf.gz.tbi"
        File parent2_vcf = "~{meta_info.members[2]}.vcf.gz"
        File parent2_vcf_index = "~{meta_info.members[2]}.vcf.gz.tbi"
        File child_gvcf = "~{meta_info.members[0]}.g.vcf.gz"
        File child_gvcf_index = "~{meta_info.members[0]}.g.vcf.gz.tbi"
        File parent1_gvcf = "~{meta_info.members[1]}.g.vcf.gz"
        File parent1_gvcf_index = "~{meta_info.members[1]}.g.vcf.gz.tbi"
        File parent2_gvcf = "~{meta_info.members[2]}.g.vcf.gz"
        File parent2_gvcf_index = "~{meta_info.members[2]}.g.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/deeptrio:1.10.0"
    }
}