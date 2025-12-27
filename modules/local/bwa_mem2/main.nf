/*
 * BWA MEM 比对模块 (低内存版本)
 */
process BWA_MEM {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    path fasta
    path "${fasta}.fai"
    path "${fasta}.amb"
    path "${fasta}.ann"
    path "${fasta}.bwt"
    path "${fasta}.pac"
    path "${fasta}.sa"
    val  output_format

    output:
    tuple val(meta), path("*.{cram,bam}"), path("*.{crai,bai}"), emit: alignment
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = meta.read_group ?: "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA"
    def format = output_format ?: 'cram'
    def output_ext = format == 'cram' ? 'cram' : 'bam'
    def sort_opts = format == 'cram' ? "--reference ${fasta} -O cram" : "-O bam"
    
    """
    bwa mem \\
        -t ${task.cpus} \\
        -R "${read_group}" \\
        ${args} \\
        ${fasta} \\
        ${reads[0]} \\
        ${reads[1]} \\
        | samtools sort -@ ${task.cpus} ${sort_opts} -o ${prefix}.${output_ext} -

    samtools index -@ ${task.cpus} ${prefix}.${output_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -oP 'Version: \\K[0-9.]+')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

/*
 * BWA-MEM2 比对模块 (高内存版本)
 */
process BWA_MEM2 {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    path fasta
    path "${fasta}.fai"
    path "${fasta}.0123"
    path "${fasta}.amb"
    path "${fasta}.ann"
    path "${fasta}.bwt.2bit.64"
    path "${fasta}.pac"
    val  output_format

    output:
    tuple val(meta), path("*.{cram,bam}"), path("*.{crai,bai}"), emit: alignment
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = meta.read_group ?: "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA"
    def format = output_format ?: 'cram'
    def output_ext = format == 'cram' ? 'cram' : 'bam'
    def sort_opts = format == 'cram' ? "--reference ${fasta} -O cram" : "-O bam"
    
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${read_group}" \\
        ${args} \\
        ${fasta} \\
        ${reads[0]} \\
        ${reads[1]} \\
        | samtools sort -@ ${task.cpus} ${sort_opts} -o ${prefix}.${output_ext} -

    samtools index -@ ${task.cpus} ${prefix}.${output_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1)
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
