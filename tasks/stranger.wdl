version 1.2

task Stranger {
    input {
        String prefix
        File vcf
        String assembly
    }

    command <<<
        if [ "~{assembly}" == "GRCh38" ]; then
            fix_assembly="grch38"
        elif [ "~{assembly}" == "GRCh37" ]; then
            fix_assembly="grch37"
        else
            echo "Unsupported assembly: ~{assembly}" >&2
            exit 1
        fi

        stranger ~{vcf} -f /stranger/resources/variant_catalog_${fix_assembly}.json > ~{prefix}.str.anno.vcf
    >>>

    output {
        File anno_vcf = "~{prefix}.str.anno.vcf"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/stranger:v0.10.0"
    }
}
