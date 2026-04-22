version 1.2

task CNVKitFix {
    input {
        String prefix
        File target_coverage
        File antitarget_coverage
        File reference
        String segmentation_method = "hmm-germline"
    }

    command <<<
        cnvkit.py fix ~{target_coverage} ~{antitarget_coverage} ~{reference} -o ~{prefix}.cnvkit.cnr
        cnvkit.py call ~{prefix}.cnvkit.cnr -o ~{prefix}.cnvkit.cns
        cnvkit.py segment ~{prefix}.cnvkit.cnr -o ~{prefix}.cnvkit.seg.cnr -m ~{segmentation_method}
        cnvkit.py call ~{prefix}.cnvkit.seg.cnr -o ~{prefix}.cnvkit.seg.cns
    >>>

    output {
        File cnvkit_cnr = "~{prefix}.cnvkit.cnr"
        File cnvkit_cns = "~{prefix}.cnvkit.cns"
        File cnvkit_seg_cnr = "~{prefix}.cnvkit.seg.cnr"
        File cnvkit_seg_cns = "~{prefix}.cnvkit.seg.cns"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/cnvkit:0.9.13.1"
    }
}

task CNVKitAntitarget {
    input {
        String prefix
        String fasta
        String assembly
        File target_bed
        Directory ref_dir
    }

    command <<<
        if [ "~{assembly}" == "GRCh38" ]; then
            excludeble_bed="/app/cnvkit/data/hg38_excludeble.bed"
        elif [ "~{assembly}" == "GRCh37" ]; then
            excludeble_bed="/app/cnvkit/data/hg19_excludeble.bed"
        else
            echo "Unsupported assembly: ~{assembly}" >&2
            exit 1
        fi

        cnvkit.py access ~{fasta} -x ${excludeble_bed} -o access.bed
        cnvkit.py antitarget ~{target_bed} -g access.bed -o ~{prefix}.antitarget.bed
    >>>

    output {
        File antitarget_bed = "~{prefix}.antitarget.bed"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/cnvkit:0.9.13.1"
    }
}

task CNVKitCoverage {
    input {
        String prefix
        File target_bed
        File antitarget_bed
        File bam
        File bai
        Int threads
    }

    command <<<
        cnvkit.py coverage ~{bam} ~{target_bed} -o ~{prefix}.targetcoverage.cnn -p ~{threads}
        cnvkit.py coverage ~{bam} ~{antitarget_bed} -o ~{prefix}.antitargetcoverage.cnn -p ~{threads}
    >>>

    output {
        File target_coverage = "~{prefix}.targetcoverage.cnn"
        File antitarget_coverage = "~{prefix}.antitargetcoverage.cnn"
    }

    runtime {
        cpu: threads
        docker: "docker.schema-bio.com/schemabio/cnvkit:0.9.13.1"
    }
}

task CNVKitReference {
    input {
        String prefix
        String fasta
        Array[File] target_coverages
        Array[File] antitarget_coverages
        Directory ref_dir
    }

    command <<<
        cnvkit.py reference -f ~{fasta} -o ~{prefix}.cnvkit.cnn \
            ~{sep=" " target_coverages} \
            ~{sep=" " antitarget_coverages}
    >>>

    output {
        File reference = "~{prefix}.cnvkit.ref.cnn"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/cnvkit:0.9.13.1"
    }
}

