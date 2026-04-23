version 1.2

task FixBed {
    input {
        File bed
    }

    command <<<
        python3 /opt/schema-germline/scripts/process_bed.py -i ~{bed} -o fixed.bed
    >>>

    output {
        File fixed_bed = "fixed.bed"
    }
}

task CreateMitoBed {
    input {
        String prefix
    }

    command <<<
        echo -e "MT\t1\t16569" > ~{prefix}.mito.bed
    >>>

    output {
        File mito_bed = "~{prefix}.mito.bed"
    }
}

task TargetBed {
    input {
        String prefix
        File bed
        String assembly
    }

    command <<<
        if [ "~{assembly}" == "GRCh38" ]; then
            excludeble_bed="/opt/schema-germline/assets/Gencode.GRCh38.cnvkit.target.bed"
        elif [ "~{assembly}" == "GRCh37" ]; then
            excludeble_bed="/opt/schema-germline/assets/Gencode.GRCh37.cnvkit.target.bed"
        else
            echo "Unsupported assembly: ~{assembly}" >&2
            exit 1
        fi

        bedtools intersect -b ~{bed} -a ${excludeble_bed} -wa -u > ~{prefix}.target.bed
    >>>

    output {
        File target_bed = "~{prefix}.target.bed"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/germline:v0.0.5"
    }
}

task FingerPrint {
    input {
        String prefix
        String fasta
        File bam
        File bai
        String assembly
        Directory ref_dir
        Int threads
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

        python3 /opt/schema-germline/scripts/sample_fingerprint.py -b ~{bam} -f ~{ref_dir}/~{fasta} \
            -s /opt/schema-germline/assets/pengelly_snp.txt -a ${fix_assembly} -t ~{threads} \
            --format json -o ~{prefix}.fingerprint.json
    >>>

    output {
        File fingerprint_json = "~{prefix}.fingerprint.json"
    }

    runtime {
        cpu: threads
        docker: "docker.schema-bio.com/schemabio/germline:v0.0.5"
    }
}
