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

        if [ "~{assembly}" == "GRCh38" ]; then
            fix_assembly="grch38"
        elif [ "~{assembly}" == "GRCh37" ]; then
            fix_assembly="grch37"
        else
            echo "Unsupported assembly: ~{assembly}" >&2
            exit 1
        fi

        ExpansionHunter \
            --reads ~{bam} \
            --reference ~{fasta} \
            --variant-catalog /app/ExpansionHunter-v5.0.0-linux_x86_64/variant_catalog/~{fix_assembly}/variant_catalog.json \
            --output-prefix ~{prefix} \
            -n ~{threads} \
            --sex ${sex}    
    >>>

    output {
        File str_vcf = "~{prefix}.str.vcf"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/expansionhunter:5.0.0"
    }
}