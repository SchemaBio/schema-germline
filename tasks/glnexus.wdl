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
# Preserve the proband caller's sites. Cohort discovery thresholds must not
# remove singleton candidates before interpretation. Parental reference blocks
# supply actual reference calls; absent coverage remains a missing genotype.
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
    }

    command <<<
        set -euo pipefail
        parents=(~{sep=" " parent_gvcfs})
        indices=(~{sep=" " parent_gvcf_tbis})
        members=(~{sep=" " sample_names})
        if [ "${#parents[@]}" -ne 2 ] || [ "${#indices[@]}" -ne 2 ] || [ "${#members[@]}" -ne 3 ]; then
            echo "ProbandFamilyVCF requires a proband and two parents" >&2
            exit 1
        fi
        reference="~{ref_dir}/~{fasta}"
        bcftools norm -f "$reference" -O z -o proband.normalized.vcf.gz "~{proband_vcf}"
        bcftools index -t proband.normalized.vcf.gz
        bcftools norm -m -any -O z -o proband.split.vcf.gz proband.normalized.vcf.gz
        bcftools index -t proband.split.vcf.gz
        for i in 0 1; do
            bcftools norm -f "$reference" -O z -o "parent.${i}.normalized.vcf.gz" "${parents[$i]}"
            bcftools index -t "parent.${i}.normalized.vcf.gz"
        done

        # Do not use --missing-to-ref. --gvcf reads actual parental reference
        # blocks, retaining no-calls when a parent's coverage is unavailable.
        bcftools merge --gvcf "$reference" -m both --threads ~{threads} -O z \
            -o family.merged.vcf.gz proband.normalized.vcf.gz \
            parent.0.normalized.vcf.gz parent.1.normalized.vcf.gz
        # Put biological identities in the declared order before indexing GT[0].
        bcftools view -s "~{sep="," sample_names}" -O z -o family.ordered.vcf.gz family.merged.vcf.gz
        expression='GT[0]="alt" && ALT!="<NON_REF>" && ALT!="<*>" && ALT!="*"'
        bcftools norm -m -any -O u family.ordered.vcf.gz \
            | bcftools view -i "$expression" -O z -o family.anchored.vcf.gz
        # Cohort QUAL/FILTER can reflect a parent's call. Restore the original
        # proband site quality, including missing values, before interpretation.
        bcftools annotate -x QUAL,FILTER -O z -o family.without-site-quality.vcf.gz family.anchored.vcf.gz
        bcftools index -t family.without-site-quality.vcf.gz
        bcftools annotate -a proband.split.vcf.gz -c QUAL,FILTER -O z \
            -o "~{prefix}.proband-family.vcf.gz" family.without-site-quality.vcf.gz
        bcftools index -t "~{prefix}.proband-family.vcf.gz"

        # Fail closed if merging lost or added a proband allele or changed its
        # genotype. Compare normalized alleles, not just chromosome/position.
        bcftools query -i "$expression" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
            proband.split.vcf.gz | LC_ALL=C sort > expected.proband.tsv
        bcftools query -s "${members[0]}" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
            "~{prefix}.proband-family.vcf.gz" | LC_ALL=C sort > actual.proband.tsv
        diff -u expected.proband.tsv actual.proband.tsv
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
