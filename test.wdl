version 1.2

workflow Test {
    input {
        String prefix
        File vcf
        String fasta
        String assembly
        Directory ref_dir
        Directory cache_dir
        Directory schema_bundle
    }

    # 参数调整
    String ref_fasta_name = basename(fasta)
    String clinvar_version = '20260415'

    # 多线程的VEP加速处理
    call GERMLINE.SplitVcf as SplitVcf {
        input:
            vcf = vcf
    }

    scatter (i in range(length(SplitVcf.split_vcfs))) {
        call VEP.VEP as VEP_Parallel {
            input:
                prefix = "~{prefix}.part~{i}",
                vcf = SplitVcf.split_vcfs[i],
                cache_dir = cache_dir,
                schema_bundle = schema_bundle,
                threads = 2,
                assembly = assembly,
                fasta = ref_fasta_name,
                clinvar_version = clinvar_version,
                ref_dir = ref_dir
        }
    }

    String merged_vcf_prefix = "~{prefix}.vep.merged"
    call GERMLINE.UniversalMergeVcfs as UniversalMergeVcfs {
        input:
            prefix = merged_vcf_prefix,
            vcfs = VEP_Parallel.out_vcf,
            threads = 8
    }
}