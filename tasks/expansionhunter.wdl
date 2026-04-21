version 1.2

task ExpansionHunter {
    input {
        String prefix
        File bam
        File bai
        String fasta
        File sry_file
        Int sex_cutoff
        Int threads
        String assembly
        Directory ref_dir
    }

    Int memory_gb = threads * 2

    command <<<
        sex="female"
        number=$(cat ~{sry_file})
        if [ $number -gt ~{sex_cutoff} ]; then
            sex="male"
        else
            sex="female"
        fi

        ExpansionHunter \
            --reads ~{bam} \
            --reference ~{fasta} \
            --variant-catalog /variant_catalog/~{assembly}/variant_catalog.json \
            --output-prefix ~{prefix} \
            -n ~{threads} \
            --sex ${sex}    
    >>>

    output {
        File str_vcf = "~{prefix}.vcf"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/expansionhunter:5.0.0"
    }
}