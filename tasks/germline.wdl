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
    }
}

# 质量统计
task QCReport {
    input {
        String prefix
        File fastp_stats
        File xamdst_json
        File mt_xamdst_json
        File fingerprint_result
        File gatk_metric
        File hs_metric
        File sry_result
        File sambamba_stats
        Int sry_cutoff
    }

    command <<<
        python /opt/schema-germline/scripts/generate_qc_report.py \
            --sample ~{prefix} \
            --output ~{prefix}.qc.json \
            --fastp ~{fastp_stats} \
            --xamdst ~{xamdst_json} \
            --mt-xamdst ~{mt_xamdst_json} \
            --fingerprint ~{fingerprint_result} \
            --metrics ~{gatk_metric} \
            --hs ~{hs_metric} \
            --sry ~{sry_result} \
            --sry-cutoff ~{sry_cutoff} \
            --sambamba-stats ~{sambamba_stats}
    >>>

    output {
        File qc_result = "~{prefix}.qc.json"
    }
}

# 转座子结果整理
task MEIReport {
    input {
        String prefix
        File mei_vcf
    }

    command <<<
        python /opt/schema-germline/scripts/mei_report.py \
            -i ~{mei_vcf} \
            -o ~{prefix}.mei.txt \
            -t /opt/schema-germline/assets/transcripts.json
    >>>

    output {
        File mei_result = "~{prefix}.mei.txt"
    }
}

# 动态突变结果整理
task STRReport {
    input {
        String prefix
        File str_vcf
    }

    command <<<
        python /opt/schema-germline/scripts/str_report.py \
            -i ~{str_vcf} \
            -o ~{prefix}.str.txt
    >>>

    output {
        File str_result = "~{prefix}.str.txt"
    }
}

# ROH突变结果整理
## GenCC版本: 20260425
task ROHReport {
    input {
        String prefix
        File automap_report
        String assembly
    }

    command <<<
        if [ "~{assembly}" == "GRCh38" ]; then
            bed=/opt/schema-germline/assets/Gencode.GRCh38.cnvkit.target.bed
        elif [ "~{assembly}" == "GRCh37" ]; then
            bed=/opt/schema-germline/assets/Gencode.GRCh37.cnvkit.target.bed
        else
            echo "Unsupported assembly: ~{assembly}" >&2
            exit 1
        fi

        python /opt/schema-germline/scripts/roh_report.py \
            -i ~{automap_report} \
            -o ~{prefix}.roh.anno.txt \
            -g /opt/schema-germline/assets/gencc-submissions.xlsx \
            -b ${bed}
    >>>

    output {
        File roh_result = "~{prefix}.roh.anno.txt"
    }
}

task SNPInDelReport {
    input {
        String prefix
        File vep_vcf
        File sry_file
        Int sex_cutoff
        String? sample_names  # Optional: comma-separated sample names (e.g., "proband,father,mother")
    }

    command <<<
        # Determine sex based on SRY reads count
        sex="female"
        number=$(cat ~{sry_file})
        if [ $number -gt ~{sex_cutoff} ]; then
            sex="male"
        else
            sex="female"
        fi

        # Set sample names (default to prefix if not provided)
        if [ -z "~{sample_names}" ]; then
            sample_names_arg="-n ~{prefix}"
        else
            sample_names_arg="-n ~{sample_names}"
        fi

        python /opt/schema-germline/scripts/vep_report.py \
            -i ~{vep_vcf} \
            -o ~{prefix}.snv_indel.txt \
            --gencc /opt/schema-germline/assets/gencc-submissions.xlsx \
            -t /opt/schema-germline/assets/transcripts.json \
            --sex ${sex} \
            ${sample_names_arg}
    >>>

    output {
        File snp_indel_result = "~{prefix}.snv_indel.txt"
    }

} 

task MTReport {
    input {
        String prefix
        File mt_vep_vcf
    }

    command <<<
        python /opt/schema-germline/scripts/mt_report.py \
            -i ~{vep_mt_vep_vcfvcf} \
            -o ~{prefix}.mt_report.txt \
            -m /opt/schema-germline/assets/mitophen.json
    >>>

    output {
        File mt_result = "~{prefix}.mt_report.txt"
    }

} 
