version 1.2

task VEP {
    input {
        String prefix
        File vcf
        Directory cache_dir
        Directory schema_bundle
        Int threads
        String assembly
        String fasta
        String clinvar_version = '20260415'
        Directory ref_dir
    }

    Int memory_gb = threads * 2

    command <<<
        cache_str="Uploaded_variation,Location,REF_ALLELE,Allele,Consequence,IMPACT,DOMAINS,Feature,DISTANCE,EXON,INTRON,SYMBOL,STRAND,HGNC_ID,HGVSc,HGVSp,HGVSg,MAX_AF,Protein_position,Amino_acids,Codons,PUBMED,Existing_variation"
        custom_str="cytoBand,CLNSIG,CLNDN,CLNSTAR,"
        gnomad_str="GnomAD_AC_joint,GnomAD_AN_joint,GnomAD_AF_joint,GnomAD_AF_joint_eas,GnomAD_nhomalt_joint_XX,GnomAD_nhomalt_joint_XY,"
        pangolin_str="Pangolin_gain_score,Pangolin_loss_score,"
        evo_str="EVOScore2_EVOScore,"
        am_str="AlphaMissense_AM,AlphaMissense_AMC,"
        self_plugin_str="FlankingSequence,MissenseZscore"

        db_prefix='hg19'
        if [ "~{assembly}" = "GRCh38" ]; then
            db_prefix='hg38'
        fi

        vep \
            --offline --cache \
            --dir_cache ~{cache_dir} --merged \
            --dir_plugins ~{cache_dir}/Plugin \
            --force_overwrite --fork ~{threads} \
            -i ~{vcf} -o ~{prefix}.vep.vcf \
            --format vcf --vcf \
            --fa ~{ref_dir}/~{fasta} \
            --shift_3prime 1 --assembly ~{assembly} --no_escape --check_existing -exclude_predicted --uploaded_allele --show_ref_allele --numbers --domains \
            --total_length --hgvs --hgvsg --symbol --ccds --uniprot --max_af --pubmed \
            --transcript_filter "stable_id match N[MR]_" \
            --plugin AnnotateClinVar,clinvar_file=~{schema_bundle}/${db_prefix}_clinvar_~{clinvar_version}.vcf.gz,fields=CLNSIG,CLNDN,CLNSTAR \
            --custom file=~{schema_bundle}/${db_prefix}_cytoBand.bed.gz,short_name=cytoBand,format=bed,type=overlap,coords=0 \
            --custom file=~{schema_bundle}/${db_prefix}_gnomad.v4.1.filtered.vcf.gz,short_name=GnomAD,format=vcf,type=exact,coords=0,fields=AC_joint%AN_joint%AF_joint%AF_joint_eas%nhomalt_joint_XX%nhomalt_joint_XY \
            --custom file=~{schema_bundle}/${db_prefix}_pangolin.vcf.gz,short_name=Pangolin,format=vcf,type=exact,coords=0,fields=gain_score%loss_score \
            --custom file=~{schema_bundle}/${db_prefix}_EVOScore2.vcf.gz,short_name=EVOScore2,format=vcf,type=exact,coords=0,fields=EVOScore \
            --custom file=~{schema_bundle}/${db_prefix}_AlphaMissense.v3.vcf.gz,short_name=AlphaMissense,format=vcf,type=exact,coords=0,fields=AM%AMC \
            --plugin FlankingSequence,10 \
            --plugin MissenseZscoreTranscript,~{schema_bundle}/missenseByTranscript.hg38.v4.1.bed \
            --fields "${cache_str},${custom_str},${gnomad_str},${pangolin_str},${evo_str},${am_str},${self_plugin_str}"
    >>>

    output {
        File out_vcf = "~{prefix}.vep.vcf"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/vep:115.2"
    }
}