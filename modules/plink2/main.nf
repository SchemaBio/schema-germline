/*
 * PLINK2 模块
 * 
 * 功能：遗传数据分析工具集
 * 工具：PLINK 2.0
 * 
 * 说明：用于基因型数据格式转换、质控、关联分析、亲缘关系推断等
 */

/*
 * VCF 转 PLINK 格式
 */
process PLINK2_VCF_TO_PLINK {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)
    path  fasta_fai  // 可选，用于设置染色体顺序

    output:
    tuple val(meta), path("*.pgen"), path("*.pvar"), path("*.psam"), emit: plink2
    tuple val(meta), path("*.bed"), path("*.bim"), path("*.fam")   , emit: plink1, optional: true
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def plink1_cmd = params.plink2_output_plink1 ?: false ? "--make-bed" : ''
    """
    plink2 \\
        --vcf ${vcf} \\
        --make-pgen \\
        ${plink1_cmd} \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version | head -n1 | sed 's/PLINK v//' | cut -d' ' -f1)
    END_VERSIONS
    """
}

/*
 * 样本/变异质控过滤
 */
process PLINK2_QC {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)

    output:
    tuple val(meta), path("*.qc.pgen"), path("*.qc.pvar"), path("*.qc.psam"), emit: plink2
    tuple val(meta), path("*.removed_samples.txt")                          , emit: removed_samples, optional: true
    tuple val(meta), path("*.removed_variants.txt")                         , emit: removed_variants, optional: true
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maf = params.plink2_maf ?: 0.01
    def geno = params.plink2_geno ?: 0.05
    def mind = params.plink2_mind ?: 0.05
    def hwe = params.plink2_hwe ?: 1e-6
    """
    plink2 \\
        --pfile ${pgen.baseName} \\
        --maf ${maf} \\
        --geno ${geno} \\
        --mind ${mind} \\
        --hwe ${hwe} \\
        --make-pgen \\
        --out ${prefix}.qc \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version | head -n1 | sed 's/PLINK v//' | cut -d' ' -f1)
    END_VERSIONS
    """
}

/*
 * 亲缘关系推断 (KING-robust)
 */
process PLINK2_RELATEDNESS {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)

    output:
    tuple val(meta), path("*.king"), path("*.king.id"), emit: king
    tuple val(meta), path("*.kin0")                   , emit: kin0, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink2 \\
        --pfile ${pgen.baseName} \\
        --make-king square \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version | head -n1 | sed 's/PLINK v//' | cut -d' ' -f1)
    END_VERSIONS
    """
}

/*
 * PCA 主成分分析
 */
process PLINK2_PCA {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)

    output:
    tuple val(meta), path("*.eigenvec"), emit: eigenvec
    tuple val(meta), path("*.eigenval"), emit: eigenval
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def n_pcs = params.plink2_pca_components ?: 10
    """
    plink2 \\
        --pfile ${pgen.baseName} \\
        --pca ${n_pcs} \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version | head -n1 | sed 's/PLINK v//' | cut -d' ' -f1)
    END_VERSIONS
    """
}

/*
 * 性别检查 (基于 X 染色体杂合率)
 */
process PLINK2_SEX_CHECK {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)

    output:
    tuple val(meta), path("*.sexcheck"), emit: sexcheck
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink2 \\
        --pfile ${pgen.baseName} \\
        --check-sex \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version | head -n1 | sed 's/PLINK v//' | cut -d' ' -f1)
    END_VERSIONS
    """
}

/*
 * 样本统计 (缺失率、杂合率)
 */
process PLINK2_SAMPLE_STATS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)

    output:
    tuple val(meta), path("*.smiss"), emit: missing
    tuple val(meta), path("*.het")  , emit: het
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink2 \\
        --pfile ${pgen.baseName} \\
        --missing sample-only \\
        --het \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version | head -n1 | sed 's/PLINK v//' | cut -d' ' -f1)
    END_VERSIONS
    """
}

/*
 * 频率统计
 */
process PLINK2_FREQ {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)

    output:
    tuple val(meta), path("*.afreq"), emit: freq
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink2 \\
        --pfile ${pgen.baseName} \\
        --freq \\
        --out ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version | head -n1 | sed 's/PLINK v//' | cut -d' ' -f1)
    END_VERSIONS
    """
}
