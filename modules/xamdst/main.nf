// XAMDST 覆盖度检查

process XAMDST {
    tag "XAMDST on $sample_id"
    label 'process_low'
    label 'mapping'
    publishDir "${params.output}/01.QC", mode: 'copy'

    input:
        tuple val(sample_id), path(alignment)
        path(bed)

    output:
        tuple val(sample_id), path("${sample_id}.xamdst.json"), emit: coverage_report

    script:
    """
    xamdst -p ${bed} -i ${alignment} -o output
    cp output/coverage.report.json ${sample_id}.xamdst.json
    """
}
