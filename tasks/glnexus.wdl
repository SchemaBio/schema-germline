version 1.2

task GLNexus {
    input {
        String prefix
        Array[File] gvcf
        Array[File] gvcf_tbi
        Int threads
    }

    Int memory_gb = threads * 2

    command <<<
        mkdir -p inputs_dir
        GVCFS=(~{sep=" " gvcf})
        TBIS=(~{sep=" " gvcf_tbi})
        for i in "${!GVCFS[@]}"; do
            ln -s "${GVCFS[$i]}" inputs_dir/
            ln -s "${TBIS[$i]}" inputs_dir/
        done

        # 注意：这里的 inputs_dir/*.g.vcf.gz 假设你的输入文件以后缀 .g.vcf.gz 结尾
        glnexus_cli --config DeepVariantWES \
            --threads ~{threads} \
            inputs_dir/*.g.vcf.gz > ~{prefix}.glnexus.bcf

        echo "Converting BCF to VCF.gz..."
        bcftools view \
            --threads ~{threads} \
            -O z \
            -o ~{prefix}.glnexus.vcf.gz \
            ~{prefix}.glnexus.bcf

        echo "Indexing VCF..."
        bcftools index -t ~{prefix}.glnexus.vcf.gz
    >>>

    output {
        File out_vcf = "~{prefix}.glnexus.vcf.gz"
        File out_tbi = "~{prefix}.glnexus.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "schemabio/glnexus:v1.4.1"
    }
}

# Proband calls define the interpreted sites; parental evidence is restricted
# to overlapping normalized loci before gVCF merging.
task ProbandFamilyVCF {
    input {
        String prefix
        File proband_vcf
        File proband_vcf_tbi
        Array[File] parent_gvcfs
        Array[File] parent_gvcf_tbis
        Array[String] sample_names
        String fasta
        Directory ref_dir
        Int threads = 4
        Int max_seconds = 120
    }

    command <<<
        set -euo pipefail
        # Bound the entire operation, including blocked I/O and all pipe stages.
        if [ ~{max_seconds} -le 0 ]; then echo "max_seconds must be positive" >&2; exit 2; fi
        timeout --signal=TERM --kill-after=10s ~{max_seconds} bash -euo pipefail <<'PROBAND_FAMILY' || { code=$?; printf '[ProbandFamilyVCF] failed exit_code=%s\n' "$code" >&2; exit "$code"; }
        log_phase() { printf '[ProbandFamilyVCF] phase=%s elapsed_seconds=%s\n' "$1" "$SECONDS" >&2; }
        parents=(~{sep=" " parent_gvcfs})
        indices=(~{sep=" " parent_gvcf_tbis})
        members=(~{sep=" " sample_names})
        if [ "${#parents[@]}" -ne 2 ] || [ "${#indices[@]}" -ne 2 ] || [ "${#members[@]}" -ne 3 ]; then
            echo "ProbandFamilyVCF requires a proband and two parents" >&2
            exit 1
        fi
        reference="~{ref_dir}/~{fasta}"
        expression='GT[0]="alt" && ALT!="<NON_REF>" && ALT!="<*>" && ALT!="*"'
        log_phase normalize-proband
        bcftools norm -f "$reference" -O z -o proband.normalized.vcf.gz "~{proband_vcf}"
        bcftools index -t proband.normalized.vcf.gz
        bcftools norm -m -any -O z -o proband.split.vcf.gz proband.normalized.vcf.gz
        bcftools index -t proband.split.vcf.gz

        # Merge overlapping BED intervals in the existing VCF contig order.
        # Include the complete REF span so deletions and gVCF END blocks overlap.
        bcftools query -i "$expression" -f '%CHROM\t%POS\t%REF\n' proband.split.vcf.gz \
            | awk 'BEGIN {OFS="\t"} {
                a=$2-1; b=a+length($3);
                if (have && $1==chrom && a<=stop) {if (b>stop) stop=b}
                else {if (have) print chrom,start,stop; chrom=$1; start=a; stop=b; have=1}
            } END {if (have) print chrom,start,stop}' > proband.regions.bed

        for i in 0 1; do
            log_phase "parent-evidence-$i"
            # One sequential pass. Indexed -R over tens of thousands of sparse
            # loci repeatedly decompresses BGZF blocks in older htslib versions.
            # Normalize before selection to preserve equivalent indel evidence.
            bcftools norm -f "$reference" -O v "${parents[$i]}" \
                | awk 'BEGIN {
                    FS=OFS="\t";
                    while ((getline line < "proband.regions.bed")>0) {
                        split(line,region,"\t"); c=region[1]; n[c]++;
                        starts[c,n[c]]=region[2]+0; stops[c,n[c]]=region[3]+0;
                    }
                    close("proband.regions.bed");
                }
                /^#/ {print; next}
                {
                    c=$1; a=$2-1; b=a+length($4);
                    if (seen[c] && a<last[c]) {print "Parental evidence is not coordinate sorted" > "/dev/stderr"; exit 2}
                    seen[c]=1; last[c]=a;
                    if (match($8,/(^|;)END=[0-9]+/)) {
                        value=substr($8,RSTART,RLENGTH); sub(/^;?END=/,"",value); b=value+0;
                    }
                    j=cursor[c]; if (!j) j=1;
                    while (j<=n[c] && stops[c,j]<=a) j++;
                    cursor[c]=j;
                    if (j<=n[c] && starts[c,j]<b) {print; kept++}
                }
                END {printf "[ProbandFamilyVCF] parental_records_selected=%d\n",kept > "/dev/stderr"}' \
                | bcftools view -O z -o "parent.${i}.evidence.vcf.gz"
            bcftools index -t "parent.${i}.evidence.vcf.gz"
        done
        log_phase merge-evidence
        # Never use --missing-to-ref. Reference genotypes require actual gVCF
        # coverage; missing evidence stays unknown. Stream BCF through reorder,
        # split and proband filtering instead of writing two full-family VCFs.
        bcftools merge --gvcf "$reference" -m both --threads ~{threads} -O u \
            proband.normalized.vcf.gz parent.0.evidence.vcf.gz parent.1.evidence.vcf.gz \
            | bcftools view -s "~{sep="," sample_names}" -O u \
            | bcftools norm -m -any -O u \
            | bcftools view -i "$expression" -O z -o family.anchored.vcf.gz

        log_phase restore-proband-quality
        bcftools annotate -x QUAL,FILTER -O z -o family.without-site-quality.vcf.gz family.anchored.vcf.gz
        bcftools index -t family.without-site-quality.vcf.gz
        bcftools annotate -a proband.split.vcf.gz -c QUAL,FILTER -O z \
            -o "~{prefix}.proband-family.vcf.gz" family.without-site-quality.vcf.gz
        bcftools index -t "~{prefix}.proband-family.vcf.gz"

        log_phase verify-proband
        bcftools query -i "$expression" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
            proband.split.vcf.gz | LC_ALL=C sort > expected.proband.tsv
        bcftools query -s "${members[0]}" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
            "~{prefix}.proband-family.vcf.gz" | LC_ALL=C sort > actual.proband.tsv
        if ! cmp -s expected.proband.tsv actual.proband.tsv; then
            echo "ProbandFamilyVCF allele/genotype preservation check failed" >&2
            exit 1
        fi
        log_phase complete
        PROBAND_FAMILY
    >>>

    output {
        File out_vcf = "~{prefix}.proband-family.vcf.gz"
        File out_tbi = "~{prefix}.proband-family.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: "8G"
        docker: "schemabio/glnexus:v1.4.1"
    }
}
