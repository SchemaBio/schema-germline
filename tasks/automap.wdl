version 1.2

task AutoMap {
    input {
        String prefix
        File vcf
        String assembly
        Int threads
    }

    Int memory_gb = threads * 2

    command <<<
        if [ "~{assembly}" == "GRCh38" ]; then
            fix_assembly="hg38"
        elif [ "~{assembly}" == "GRCh37" ]; then
            fix_assembly="hg19"
        else
            echo "Unsupported assembly: ~{assembly}" >&2
            exit 1
        fi
        
        echo -e "#Chr\tBegin\tEnd\tSize(Mb)\tNb_variants\tPercentage_homozygosity\n" > ~{prefix}.ROH.txt

        bash /opt/AutoMap/AutoMap_v1.3.sh \
            --vcf ~{vcf} \
            --genome ${fix_assembly} \
            --out result \
            --id ~{prefix}
        
        if [ -f result/~{prefix}.HomRegions.tsv ]; then
            cat result/~{prefix}.HomRegions.tsv | grep -v '##' > ~{prefix}.ROH.txt
        fi
    >>>

    output {
        File combined_tsv = "~{prefix}.ROH.txt"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/automap:1.3"
        cpu: threads
        memory: "~{memory_gb}G"
    }
}