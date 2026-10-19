version 1.2

task InitializeVcf {
    input {
        String prefix
        File vcf
        String assembly
        String fasta
        Directory ref_dir
        String? sample_name
        Int threads = 4
    }

    command <<<
        set -euo pipefail
        python3 /opt/schema-germline/scripts/initialize_vcf.py \
            --prefix "~{prefix}" \
            --vcf "~{vcf}" \
            --assembly "~{assembly}" \
            --reference "~{ref_dir}/~{fasta}" \
            --sample-name "~{default='' sample_name}" \
            --threads ~{threads}
    >>>

    output {
        File normalized_vcf = "~{prefix}.normalized.vcf.gz"
        File normalized_vcf_tbi = "~{prefix}.normalized.vcf.gz.tbi"
        File selected_vcf = "~{prefix}.selected.vcf.gz"
        File preprocessing = "~{prefix}.preprocessing.json"
        Int variant_count = read_int("variant_count.txt")
        String selected_sample = read_string("sample_name.txt")
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/germline:v0.1.7"
        cpu: threads
        memory: "8G"
    }
}

task InferVcfSex {
    input {
        String prefix
        File vcf
        String sample_name
        Int minimum_sites = 100
        Float female_max_ratio = 0.5
        Float male_min_ratio = 0.7
    }

    command <<<
        set -euo pipefail
        python3 /opt/schema-germline/scripts/infer_vcf_sex.py \
            --prefix "~{prefix}" \
            --vcf "~{vcf}" \
            --sample-name "~{sample_name}" \
            --minimum-sites ~{minimum_sites} \
            --female-max-ratio ~{female_max_ratio} \
            --male-min-ratio ~{male_min_ratio}
    >>>

    output {
        String sex = read_string("inferred_sex.txt")
        File sex_inference = "~{prefix}.sex_inference.json"
        File sex_explanation = "~{prefix}.sex_inference.txt"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/germline:v0.1.7"
        cpu: 1
        memory: "2G"
    }
}

task IndexAnnotatedVcf {
    input {
        String prefix
        File normalized_vcf
        File vep_vcf
        Int variant_count
    }

    command <<<
        set -euo pipefail
        # VEP must not silently discard or change normalized variants or GTs.
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n' "~{normalized_vcf}" | sort > expected.tsv
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n' "~{vep_vcf}" | sort > observed.tsv
        cmp expected.tsv observed.tsv
        if [ ~{variant_count} -gt 0 ]; then
            bcftools view -h "~{vep_vcf}" > annotated.header
            grep -q '^##INFO=<ID=CSQ,' annotated.header
        fi
        bcftools sort -Oz -o "~{prefix}.vep.vcf.gz" "~{vep_vcf}"
        bcftools index --tbi --force "~{prefix}.vep.vcf.gz"
    >>>

    output {
        File annotated_vcf = "~{prefix}.vep.vcf.gz"
        File annotated_vcf_tbi = "~{prefix}.vep.vcf.gz.tbi"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/germline:v0.1.7"
        cpu: 2
        memory: "4G"
    }
}
