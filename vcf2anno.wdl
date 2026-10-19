version 1.2

import "tasks/vcf2anno.wdl" as VCF
import "tasks/vep.wdl" as VEP
import "tasks/germline.wdl" as GERMLINE

workflow VCF2Anno {
    input {
        String prefix
        File vcf
        String assembly
        String fasta
        Directory ref_dir
        Directory cache_dir
        Directory schema_bundle
        String? sample_name
        Int threads = 8
        String clinvar_version = "20260415"
        Int sex_min_sites = 100
        Float female_max_ratio = 0.5
        Float male_min_ratio = 0.7
    }

    String ref_fasta_name = basename(fasta)

    call VCF.InitializeVcf as InitializeVcf {
        input:
            prefix = prefix,
            vcf = vcf,
            assembly = assembly,
            fasta = ref_fasta_name,
            ref_dir = ref_dir,
            sample_name = sample_name,
            threads = if threads > 4 then 4 else threads
    }

    # Infer from selected original records, before multiallelic splitting.
    call VCF.InferVcfSex as InferSex {
        input:
            prefix = prefix,
            vcf = InitializeVcf.selected_vcf,
            sample_name = InitializeVcf.selected_sample,
            minimum_sites = sex_min_sites,
            female_max_ratio = female_max_ratio,
            male_min_ratio = male_min_ratio
    }

    # Fork workers replace the old chromosome split/merge path for this workflow.
    if (InitializeVcf.variant_count > 0) {
        call VEP.VEP as Annotate {
            input:
                prefix = prefix,
                vcf = InitializeVcf.normalized_vcf,
                assembly = assembly,
                fasta = ref_fasta_name,
                ref_dir = ref_dir,
                cache_dir = cache_dir,
                schema_bundle = schema_bundle,
                threads = threads,
                clinvar_version = clinvar_version
        }
    }

    call VCF.IndexAnnotatedVcf as IndexAnnotatedVcf {
        input:
            prefix = prefix,
            normalized_vcf = InitializeVcf.normalized_vcf,
            vep_vcf = select_first([Annotate.out_vcf, InitializeVcf.normalized_vcf]),
            variant_count = InitializeVcf.variant_count
    }

    call GERMLINE.SNPInDelReport as SNPInDelReport {
        input:
            prefix = prefix,
            vep_vcf = IndexAnnotatedVcf.annotated_vcf,
            sex = InferSex.sex,
            external_vcf = true,
            expected_variant_count = InitializeVcf.variant_count
    }

    output {
        File normalized_vcf = InitializeVcf.normalized_vcf
        File normalized_vcf_tbi = InitializeVcf.normalized_vcf_tbi
        File annotated_vcf = IndexAnnotatedVcf.annotated_vcf
        File annotated_vcf_tbi = IndexAnnotatedVcf.annotated_vcf_tbi
        File snp_indel = SNPInDelReport.snp_indel_result
        File preprocessing = InitializeVcf.preprocessing
        File sex_inference = InferSex.sex_inference
        File sex_explanation = InferSex.sex_explanation
        File summary = write_json(object {
            pipeline: "VCF2Anno", version: "v0.1.0", status: "Success",
            prefix: prefix, assembly: assembly, sample_name: InitializeVcf.selected_sample,
            vep_version: "116.0", clinvar_version: clinvar_version,
            vep_skipped: InitializeVcf.variant_count == 0,
            normalized_vcf: InitializeVcf.normalized_vcf,
            normalized_vcf_tbi: InitializeVcf.normalized_vcf_tbi,
            annotated_vcf: IndexAnnotatedVcf.annotated_vcf,
            annotated_vcf_tbi: IndexAnnotatedVcf.annotated_vcf_tbi,
            snp_indel: SNPInDelReport.snp_indel_result,
            preprocessing: read_json(InitializeVcf.preprocessing),
            sex_inference: read_json(InferSex.sex_inference),
            sex_explanation: InferSex.sex_explanation
        })
    }
}
