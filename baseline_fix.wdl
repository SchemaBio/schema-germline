version 1.2

import "tasks/fastp.wdl" as FASTP
import "tasks/germline.wdl" as GERMLINE
import "tasks/cnvkit.wdl" as CNVKIT
import "tasks/bwamem.wdl" as BWAMEM
import "tasks/sambamba.wdl" as SAMBAMBA

struct PipelineSummary {
    String prefix
    String status
    String pipeline
    String version
    String update_method
    File reference
    File update_report
    Int previous_number_of_samples
    Int added_number_of_samples
    Int total_number_of_samples
    Array[String] added_sample_names
}

workflow CNVBaselineFix {
    input {
        String prefix
        File bed
        String fasta
        String assembly
        File existing_reference
        Int existing_sample_count
        Array[File] read_1
        Array[File] read_2
        Directory ref_dir
    }

    String ref_fasta_name = basename(fasta)
    Int total_threads = 32

    call GERMLINE.ValidateBaselineFixInputs as ValidateInputs {
        input:
            read_1_count = length(read_1),
            read_2_count = length(read_2),
            existing_sample_count = existing_sample_count
    }

    Int added_sample_count = ValidateInputs.added_sample_count
    Int base_threads = total_threads / added_sample_count
    Int fastp_threads = if base_threads < 8 then 8 else if base_threads > 16 then 16 else base_threads
    Int bwa_threads = if base_threads < 16 then 16 else if base_threads > 32 then 32 else base_threads
    Int cnvkit_threads = if base_threads < 8 then 8 else if base_threads > 16 then 16 else base_threads

    call GERMLINE.FixBed as FixBed {
        input:
            bed = bed
    }

    call GERMLINE.TargetBed as TargetBed {
        input:
            prefix = prefix,
            bed = FixBed.fixed_bed,
            assembly = assembly
    }

    call CNVKIT.CNVKitAntitarget as CNVKitAntitarget {
        input:
            prefix = prefix,
            target_bed = TargetBed.target_bed,
            fasta = ref_fasta_name,
            assembly = assembly,
            ref_dir = ref_dir
    }

    scatter (i in range(added_sample_count)) {
        String raw_filename = basename(read_1[i])
        String clean_name = sub(raw_filename, "\\.fastq\\.gz$|\\.fq\\.gz$", "")
        String sample_prefix = "~{prefix}_new_sample~{i}"

        call FASTP.Fastp as Fastp {
            input:
                prefix = sample_prefix,
                read_1 = read_1[i],
                read_2 = read_2[i],
                threads = fastp_threads
        }

        call BWAMEM.BwaAlign as BwaAlign {
            input:
                prefix = sample_prefix,
                read_1 = Fastp.clean_read_1,
                read_2 = Fastp.clean_read_2,
                ref_dir = ref_dir,
                ref_fasta_name = ref_fasta_name,
                threads = bwa_threads
        }

        call SAMBAMBA.SambambaMarkdup as Markdup {
            input:
                prefix = sample_prefix,
                bam = BwaAlign.out_bam,
                bai = BwaAlign.out_bai,
                threads = bwa_threads
        }

        call CNVKIT.CNVKitCoverage as CNVKitCoverage {
            input:
                prefix = sample_prefix,
                target_bed = TargetBed.target_bed,
                antitarget_bed = CNVKitAntitarget.antitarget_bed,
                bam = Markdup.markdup_bam,
                bai = Markdup.markdup_bai,
                threads = cnvkit_threads
        }
    }

    call CNVKIT.CNVKitReference as NewBatchReference {
        input:
            prefix = "~{prefix}.new_batch",
            fasta = ref_fasta_name,
            target_coverages = CNVKitCoverage.target_coverage,
            antitarget_coverages = CNVKitCoverage.antitarget_coverage,
            ref_dir = ref_dir
    }

    call GERMLINE.UpdateCNVKitReference as UpdateReference {
        input:
            prefix = prefix,
            old_reference = existing_reference,
            new_reference = NewBatchReference.reference,
            old_sample_count = existing_sample_count,
            new_sample_count = added_sample_count
    }

    Int total_sample_count = existing_sample_count + added_sample_count

    output {
        File reference = UpdateReference.reference
        File update_report = UpdateReference.update_report
        File new_batch_reference = NewBatchReference.reference
        Array[File] new_target_coverages = CNVKitCoverage.target_coverage
        Array[File] new_antitarget_coverages = CNVKitCoverage.antitarget_coverage
        File summary = write_json(
            PipelineSummary {
                prefix: prefix,
                status: "Success",
                pipeline: "CNV_Baseline_Fix",
                version: "v0.0.1",
                update_method: "sample_count_weighted_reference_merge",
                reference: UpdateReference.reference,
                update_report: UpdateReference.update_report,
                previous_number_of_samples: existing_sample_count,
                added_number_of_samples: added_sample_count,
                total_number_of_samples: total_sample_count,
                added_sample_names: clean_name
            }
        )
    }
}
