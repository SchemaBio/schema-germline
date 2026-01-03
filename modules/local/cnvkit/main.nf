/*
 * CNVkit 模块集合
 *
 * 功能：基于测序深度的拷贝数变异 (CNV) 检测
 * 工具：CNVkit
 *
 * 包含：
 *   - CNVKIT_REFERENCE: 构建参考基线 (使用正常样本，内含 coverage 计算)
 *   - CNVKIT_REFERENCE_FLAT: 构建平坦基线 (无正常样本)
 *   - CNVKIT_CALL: CNV 检测 (coverage + fix + segment + call)
 *   - CNVKIT_EXPORT_VCF: 导出 VCF 格式
 *   - CNVKIT_SEX: 性别检测
 */

/*
 * CNVKIT_REFERENCE - 构建参考基线
 *
 * 功能：使用多个正常样本构建参考基线 (内含 coverage 计算)
 * 输入：多个 BAM/CRAM 文件
 * 输出：reference.cnn 基线文件
 *
 * 说明：至少需要 1 个正常样本，推荐 >= 5 个以获得稳定基线
 */
process CNVKIT_REFERENCE {
    tag "reference"
    label 'process_high'

    input:
    path alignments      // 多个 BAM/CRAM 文件
    path alignment_indexes
    path fasta
    path fasta_fai
    path targets      // 目标区域 BED
    path antitargets  // 反目标区域 BED

    output:
    path "reference.cnn", emit: reference
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ref_cmd = alignments[0].name.endsWith('.cram') ? "--fasta ${fasta}" : ''
    """
    # 为每个样本计算 coverage
    for bam in ${alignments}; do
        sample_name=\$(basename \${bam} | sed 's/\\.[^.]*\$//')

        cnvkit.py coverage \\
            \${bam} \\
            ${targets} \\
            ${ref_cmd} \\
            --processes ${task.cpus} \\
            --output \${sample_name}.targetcoverage.cnn

        cnvkit.py coverage \\
            \${bam} \\
            ${antitargets} \\
            ${ref_cmd} \\
            --processes ${task.cpus} \\
            --output \${sample_name}.antitargetcoverage.cnn
    done

    # 构建参考基线
    cnvkit.py reference \\
        *.targetcoverage.cnn \\
        *.antitargetcoverage.cnn \\
        --fasta ${fasta} \\
        --output reference.cnn \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed 's/cnvkit //')
    END_VERSIONS
    """
}

/*
 * CNVKIT_REFERENCE_FLAT - 构建平坦参考基线 (无正常样本)
 *
 * 功能：当没有正常样本时，使用目标区域文件构建平坦基线
 * 说明：平坦基线假设所有区域拷贝数为 2，检测效果不如使用正常样本
 */
process CNVKIT_REFERENCE_FLAT {
    tag "flat_reference"
    label 'process_low'

    input:
    path targets      // 目标区域 BED 文件
    path antitargets  // 反目标区域 BED 文件
    path fasta
    path fasta_fai

    output:
    path "flat_reference.cnn", emit: reference
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cnvkit.py reference \\
        --fasta ${fasta} \\
        --targets ${targets} \\
        --antitargets ${antitargets} \\
        --output flat_reference.cnn \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed 's/cnvkit //')
    END_VERSIONS
    """
}

/*
 * CNVKIT_CALL - CNV 检测
 *
 * 功能：对单个样本进行 CNV 检测 (内含 coverage 计算)
 * 步骤：coverage -> fix -> segment -> call
 * 输出：
 *   - .cnr: 修正后的拷贝数比值
 *   - .cns: 分段结果
 *   - .call.cns: CNV 调用结果
 */
process CNVKIT_CALL {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(alignment), path(alignment_index)
    path  fasta
    path  fasta_fai
    path  targets      // 目标区域 BED
    path  antitargets  // 反目标区域 BED
    path  reference    // 参考基线 .cnn 文件

    output:
    tuple val(meta), path("*.cnr")     , emit: cnr
    tuple val(meta), path("*.cns")     , emit: cns
    tuple val(meta), path("*.call.cns"), emit: call
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args_fix = task.ext.args_fix ?: ''
    def args_segment = task.ext.args_segment ?: ''
    def args_call = task.ext.args_call ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_cmd = alignment.name.endsWith('.cram') ? "--fasta ${fasta}" : ''
    def ploidy = meta.sex == 'male' ? '--male-reference' : ''
    """
    # Step 1: Coverage - 计算覆盖度
    cnvkit.py coverage \\
        ${alignment} \\
        ${targets} \\
        ${ref_cmd} \\
        --processes ${task.cpus} \\
        --output ${prefix}.targetcoverage.cnn

    cnvkit.py coverage \\
        ${alignment} \\
        ${antitargets} \\
        ${ref_cmd} \\
        --processes ${task.cpus} \\
        --output ${prefix}.antitargetcoverage.cnn

    # Step 2: Fix - 修正覆盖度偏差
    cnvkit.py fix \\
        ${prefix}.targetcoverage.cnn \\
        ${prefix}.antitargetcoverage.cnn \\
        ${reference} \\
        --output ${prefix}.cnr \\
        ${args_fix}

    # Step 3: Segment - CBS 分段
    cnvkit.py segment \\
        ${prefix}.cnr \\
        --processes ${task.cpus} \\
        --output ${prefix}.cns \\
        ${args_segment}

    # Step 4: Call - CNV 调用
    cnvkit.py call \\
        ${prefix}.cns \\
        ${ploidy} \\
        --output ${prefix}.call.cns \\
        ${args_call}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed 's/cnvkit //')
    END_VERSIONS
    """
}

/*
 * CNVKIT_EXPORT_VCF - 导出 VCF 格式
 */
process CNVKIT_EXPORT_VCF {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(call_cns)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ploidy = meta.sex == 'male' ? '--male-reference' : ''
    """
    cnvkit.py export vcf \\
        ${call_cns} \\
        ${ploidy} \\
        --output ${prefix}.cnv.vcf \\
        ${args}

    bgzip -c ${prefix}.cnv.vcf > ${prefix}.cnv.vcf.gz
    tabix -p vcf ${prefix}.cnv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed 's/cnvkit //')
    END_VERSIONS
    """
}

/*
 * CNVKIT_SEX - 性别检测
 *
 * 功能：根据 X/Y 染色体覆盖度推断样本性别
 * 输入：BAM/CRAM 文件 + targets/antitargets BED
 * 输出：性别推断结果
 */
process CNVKIT_SEX {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(alignment), path(alignment_index)
    path  fasta
    path  fasta_fai
    path  targets
    path  antitargets

    output:
    tuple val(meta), path("*.sex_check.txt"), emit: result
    tuple val(meta), env(INFERRED_SEX)      , emit: sex
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_cmd = alignment.name.endsWith('.cram') ? "--fasta ${fasta}" : ''
    def declared_sex = meta.sex ?: 'unknown'
    """
    # 计算覆盖度
    cnvkit.py coverage \\
        ${alignment} \\
        ${targets} \\
        ${ref_cmd} \\
        --processes ${task.cpus} \\
        --output ${prefix}.targetcoverage.cnn

    cnvkit.py coverage \\
        ${alignment} \\
        ${antitargets} \\
        ${ref_cmd} \\
        --processes ${task.cpus} \\
        --output ${prefix}.antitargetcoverage.cnn

    # 性别推断
    SEX_RESULT=\$(cnvkit.py sex ${prefix}.targetcoverage.cnn ${prefix}.antitargetcoverage.cnn ${args})
    INFERRED_SEX=\$(echo "\${SEX_RESULT}" | awk '{print \$2}')

    # 转换性别标识
    if [ "\${INFERRED_SEX}" == "Male" ]; then
        INFERRED_SEX="male"
    elif [ "\${INFERRED_SEX}" == "Female" ]; then
        INFERRED_SEX="female"
    fi

    # 一致性检查
    DECLARED="${declared_sex}"
    if [ "\${DECLARED}" != "unknown" ] && [ "\${DECLARED}" != "\${INFERRED_SEX}" ]; then
        STATUS="MISMATCH"
    else
        STATUS="PASS"
    fi

    # 输出结果
    cat <<-EOF > ${prefix}.sex_check.txt
    sample_id: ${meta.id}
    inferred_sex: \${INFERRED_SEX}
    declared_sex: ${declared_sex}
    status: \${STATUS}
    raw_output: \${SEX_RESULT}
    EOF

    export INFERRED_SEX

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed 's/cnvkit //')
    END_VERSIONS
    """
}