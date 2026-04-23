version 1.2

task Stranger {
    input {
        String prefix
        File vcf
        String assembly
    }

    # Convert assembly format in WDL, not shell
    String fix_assembly = if assembly == "GRCh38" then "grch38"
                         else if assembly == "GRCh37" then "grch37"
                         else "grch38"  # default fallback

    command {
        stranger ~{vcf} -f /app/stranger/stranger/resources/variant_catalog_~{fix_assembly}.json > ~{prefix}.str.anno.vcf
    }

    output {
        File anno_vcf = "~{prefix}.str.anno.vcf"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/stranger:v0.10.0.1"
    }
}
