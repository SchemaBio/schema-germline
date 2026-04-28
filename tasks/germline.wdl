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
        python3 /opt/schema-germline/scripts/generate_qc_report.py \
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
        python3 /opt/schema-germline/scripts/mei_report.py \
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
        python3 /opt/schema-germline/scripts/str_report.py \
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

        python3 /opt/schema-germline/scripts/roh_report.py \
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

        python3 /opt/schema-germline/scripts/vep_report.py \
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
        python3 /opt/schema-germline/scripts/mt_report.py \
            -i ~{mt_vep_vcf} \
            -o ~{prefix}.mt_report.txt \
            -m /opt/schema-germline/assets/mitophen.json
    >>>

    output {
        File mt_result = "~{prefix}.mt_report.txt"
    }

} 

task CNVGene {
    input {
        String prefix
        File cnv_cnr
        Float cnv_del_threshold = 1.5
        Float cnv_dup_threshold = 2.5
    }

    command <<<
        python3 /opt/schema-germline/scripts/gene_level_segment.py \
            -i ~{cnv_cnr} \
            -o ~{prefix}.cnv.gene.txt \
            --cn-del-threshold ~{cnv_del_threshold} \
            --cn-dup-threshold ~{cnv_dup_threshold}
    >>>

    output {
        File cnv_result = "~{prefix}.cnv.gene.txt"
    }

}

task CNVRegion {
    input {
        String prefix
        File cnv_cnr
        Int bin_size
        Float cnv_del_threshold = 1.5
        Float cnv_dup_threshold = 2.5
    }

    command <<<
        python3 /opt/schema-germline/scripts/merge_bins_for_cnv.py \
            -i ~{cnv_cnr} \
            -o ~{prefix}.cnv.region.txt \
            --del-threshold ~{cnv_del_threshold} \
            --dup-threshold ~{cnv_dup_threshold} \
            --bin-size ~{bin_size} \
            --keep-antitarget
    >>>

    output {
        File cnv_result = "~{prefix}.cnv.region.txt"
    }

} 

task UPD {
    input {
        String prefix
        File proband_roh_report
        File father_roh_report
        File mother_roh_report
    }

    command <<<
        python3 /opt/schema-germline/scripts/upd_detection.py \
            -p ~{proband_roh_report} \
            -f ~{father_roh_report} \
            -m ~{mother_roh_report} \
            -o ~{prefix}.upd_report.txt
    >>>

    output {
        File upd_result = "~{prefix}.upd_report.txt"
    }

} 

# 拆分任务
task SplitVcf {
    input {
        File vcf
    }
    
    command <<<
        set -ex
        # 确保输入 VCF 有索引（如果输入没给，就现场建一个）
        if [ ! -f "~{vcf}.tbi" ]; then
            bcftools index -t ~{vcf}
        fi

        for chrom in $(bcftools view -h ~{vcf} | grep "^##contig" | sed 's/.*ID=\([^,>]*\).*/\1/' | grep -E "^(chr)?([0-9]+|X|Y|M)$"); do
            bcftools view ~{vcf} --regions $chrom -O z -o ${chrom}.split.vcf.gz
            bcftools index -t ${chrom}.split.vcf.gz
        done
    >>>
    
    output {
        Array[File] split_vcfs = glob("*.split.vcf.gz")
        Array[File] split_vcf_tbis = glob("*.split.vcf.gz.tbi")
    }
}

# 合并任务
task UniversalMergeVcfs {
    input {
        String prefix
        Array[File] vcfs
        Int threads
    }

    Int memory_gb = threads * 4

    command <<<
        set -ex

        # 1. 汇聚文件并建立软链接
        # 将分散在不同目录的 VCF 和 TBI 链接到当前目录，确保 bcftools 能找到配对索引
        VCF_ARRAY=(~{sep=' ' vcfs})

        for f in "${VCF_ARRAY[@]}"; do
            fname=$(basename "$f")
            # 检查是否为压缩格式 (.vcf.gz 或 .vcf.bgz)
            if [[ "$fname" == *.vcf.gz ]] || [[ "$fname" == *.vcf.bgz ]]; then
                ln -sf "$f" "$fname"
            else
                # 未压缩的 VCF，用 bcftools 压缩
                ln -sf "$f" "$fname"
                compressed_name="${fname}.gz"
                bcftools view -O z -o "$compressed_name" "$fname"
                fname="$compressed_name"
            fi
            bcftools index -t "$fname"
            echo "$fname" >> processed_vcfs.txt
        done

        # 2. 生成基于基因组坐标的有序列表
        while read -r f; do
            INFO=$(bcftools view -H "$f" | head -n 1 | awk '{print $1"\t"$2}')
            if [ -n "$INFO" ]; then
                printf "%s\t%s\n" "$f" "$INFO" >> unsorted_list.tmp
            fi
        done < processed_vcfs.txt

        # 排序逻辑：
        sort -k2,2V -k3,3n unsorted_list.tmp | cut -f1 > final_vcf_list.txt

        # 3. 使用文件列表执行合并
        bcftools concat -f final_vcf_list.txt -a -O z -o ~{prefix}.merged.vcf.gz

        # 4. 建立最终索引
        bcftools index -t ~{prefix}.merged.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.merged.vcf.gz"
        File merged_vcf_tbi = "~{prefix}.merged.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
    }
}