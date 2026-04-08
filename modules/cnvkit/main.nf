// CNVKIT 模块
// 用途：基于测序深度的拷贝数变异 (CNV) 检测
// 工具：CNVkit
//
// 分析流程 (单样本模式，已有基线文件)：
//   1. CNVKIT_TARGET - 从用户 BED 生成 target BED
//   2. CNVKIT_ANTITARGET - 从 target BED 生成 antitarget BED
//   3. CNVKIT_COVERAGE - 计算覆盖度
//   4. CNVKIT_FIX - 修正覆盖度偏差 (需要基线文件)
//   5. CNVKIT_SEGMENT - CBS 分段
//   6. CNVKIT_CALL - CNV 调用
//   7. CNVKIT_EXPORT_VCF - 导出 VCF 格式
//
// 说明：
//   - 需要预先准备基线文件 (reference.cnn)
//   - 用户 BED 文件会自动转换为 target/antitarget 格式
//   - 支持自动性别检测和性别特异分析

// ============================================================================
// CNVKIT_TARGET - 生成目标区域 BED
// ============================================================================
process CNVKIT_TARGET {
    tag "CNVKIT_TARGET"
    label 'process_low'
    label 'cnvkit'
    publishDir "${params.output}/04.CNV/target", mode: 'copy'

    input:
        path user_bed           // 用户提供的捕获区域 BED 文件
        val annotate            // 是否添加基因注释 (boolean)
        val split_size          // 分割大小 (默认 5000)

    output:
        path "targets.bed", emit: target_bed

    script:
    def split = split_size ?: 5000
    def annotate_cmd = annotate ? '--annotate' : ''
    """
    cnvkit.py target \\
        ${user_bed} \\
        --output targets.bed \\
        --split ${split} \\
        ${annotate_cmd}
    """
}

// ============================================================================
// CNVKIT_ANTITARGET - 生成反目标区域 BED
// ============================================================================
process CNVKIT_ANTITARGET {
    tag "CNVKIT_ANTITARGET"
    label 'process_low'
    label 'cnvkit'
    publishDir "${params.output}/04.CNV/target", mode: 'copy'

    input:
        path target_bed         // 目标区域 BED
        path access_bed         // 可访问区域 BED (可选，用于限制 antitarget 范围)
        val min_target_size     // 最小目标大小 (默认 10000)

    output:
        path "antitargets.bed", emit: antitarget_bed

    script:
    def min_size = min_target_size ?: 10000
    def access_cmd = access_bed.name != 'NO_FILE' ? "--access ${access_bed}" : ''
    """
    cnvkit.py antitarget \\
        ${target_bed} \\
        --output antitargets.bed \\
        --min-target-size ${min_size} \\
        ${access_cmd}
    """
}

// ============================================================================
// CNVKIT_COVERAGE - 计算覆盖度
// ============================================================================
process CNVKIT_COVERAGE {
    tag "CNVKIT_COVERAGE on ${alignment.baseName}"
    label 'process_medium'
    label 'cnvkit'
    publishDir "${params.output}/04.CNV/coverage", mode: 'copy'

    input:
        path alignment          // BAM/CRAM 比对文件
        path alignment_index    // 比对文件索引
        path fasta              // 参考基因组 (CRAM 输入时需要)
        path target_bed         // 目标区域 BED
        path antitarget_bed     // 反目标区域 BED

    output:
        path "*.targetcoverage.cnn", emit: target_coverage
        path "*.antitargetcoverage.cnn", emit: antitarget_coverage

    script:
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam|cram)$/, '')
    def ref_cmd = alignment.name.endsWith('.cram') ? "--fasta ${fasta}" : ''
    """
    # 计算 target 区域覆盖度
    cnvkit.py coverage \\
        ${alignment} \\
        ${target_bed} \\
        ${ref_cmd} \\
        --processes ${task.cpus} \\
        --output ${sample_id}.targetcoverage.cnn

    # 计算 antitarget 区域覆盖度
    cnvkit.py coverage \\
        ${alignment} \\
        ${antitarget_bed} \\
        ${ref_cmd} \\
        --processes ${task.cpus} \\
        --output ${sample_id}.antitargetcoverage.cnn
    """
}

// ============================================================================
// CNVKIT_FIX - 修正覆盖度偏差
// ============================================================================
process CNVKIT_FIX {
    tag "CNVKIT_FIX on ${target_coverage.baseName}"
    label 'process_low'
    label 'cnvkit'
    publishDir "${params.output}/04.CNV/fix", mode: 'copy'

    input:
        path target_coverage    // 目标区域覆盖度文件 (.targetcoverage.cnn)
        path antitarget_coverage // 反目标区域覆盖度文件 (.antitargetcoverage.cnn)
        path reference          // 基线文件 (reference.cnn)
        path target_bed         // 目标区域 BED (用于 QC)

    output:
        path "*.cnr", emit: cnr

    script:
    def sample_id = target_coverage.baseName.replaceAll(/\\.targetcoverage$/, '')
    """
    cnvkit.py fix \\
        ${target_coverage} \\
        ${antitarget_coverage} \\
        ${reference} \\
        --output ${sample_id}.cnr
    """
}

// ============================================================================
// CNVKIT_SEGMENT - CBS 分段
// ============================================================================
process CNVKIT_SEGMENT {
    tag "CNVKIT_SEGMENT on ${cnr.baseName}"
    label 'process_medium'
    label 'cnvkit'
    publishDir "${params.output}/04.CNV/segment", mode: 'copy'

    input:
        path cnr                // 修正后的覆盖度文件
        val method              // 分段方法: 'cbs' (默认), 'haar', 'flasso'
        val threshold           // 分段阈值 (默认自动)

    output:
        path "*.cns", emit: cns

    script:
    def seg_method = method ?: 'cbs'
    def thresh_cmd = threshold ? "--threshold ${threshold}" : ''
    def sample_id = cnr.baseName
    """
    cnvkit.py segment \\
        ${cnr} \\
        --method ${seg_method} \\
        --processes ${task.cpus} \\
        ${thresh_cmd} \\
        --output ${sample_id}.cns
    """
}

// ============================================================================
// CNVKIT_CALL - CNV 调用
// ============================================================================
process CNVKIT_CALL {
    tag "CNVKIT_CALL on ${cns.baseName}"
    label 'process_low'
    label 'cnvkit'
    publishDir "${params.output}/04.CNV/call", mode: 'copy'

    input:
        path cns                // 分段结果文件
        val ploidy              // 倍性 (默认 2)
        val purity              // 肿瘤纯度 (可选，0-1)
        val is_male             // 样本是否为男性 (boolean)
        val drop_outliers       // 是否丢弃离群值 (boolean, 默认 true)

    output:
        path "*.call.cns", emit: call_cns

    script:
    def sample_id = cns.baseName
    def ploidy_val = ploidy ?: 2
    def purity_cmd = purity ? "--purity ${purity}" : ''
    def male_cmd = is_male ? '--male-reference' : ''
    def drop_cmd = drop_outliers ? '--drop-low-coverage' : ''
    """
    cnvkit.py call \\
        ${cns} \\
        --ploidy ${ploidy_val} \\
        ${purity_cmd} \\
        ${male_cmd} \\
        ${drop_cmd} \\
        --output ${sample_id}.call.cns
    """
}

// ============================================================================
// CNVKIT_EXPORT_VCF - 导出 VCF 格式
// ============================================================================
process CNVKIT_EXPORT_VCF {
    tag "CNVKIT_EXPORT_VCF on ${call_cns.baseName}"
    label 'process_low'
    label 'cnvkit'
    publishDir "${params.output}/04.CNV/vcf", mode: 'copy'

    input:
        path call_cns           // CNV 调用结果
        val is_male             // 样本是否为男性 (boolean)

    output:
        path "*.cnv.vcf.gz", emit: vcf
        path "*.cnv.vcf.gz.tbi", emit: vcf_tbi

    script:
    def sample_id = call_cns.baseName.replaceAll(/\\.call$/, '')
    def male_cmd = is_male ? '--male-reference' : ''
    """
    # 导出 VCF 格式
    cnvkit.py export vcf \\
        ${call_cns} \\
        ${male_cmd} \\
        --output ${sample_id}.cnv.vcf

    # 压缩并索引
    bgzip -c ${sample_id}.cnv.vcf > ${sample_id}.cnv.vcf.gz
    tabix -p vcf ${sample_id}.cnv.vcf.gz

    # 清理临时文件
    rm -f ${sample_id}.cnv.vcf
    """
}

// ============================================================================
// CNVKIT_SEX - 性别检测
// ============================================================================
process CNVKIT_SEX {
    tag "CNVKIT_SEX on ${target_coverage.baseName}"
    label 'process_low'
    label 'cnvkit'
    publishDir "${params.output}/04.CNV/sex_check", mode: 'copy'

    input:
        path target_coverage    // 目标区域覆盖度文件
        path antitarget_coverage // 反目标区域覆盖度文件

    output:
        path "*.sex.json", emit: sex_json

    script:
    def sample_id = target_coverage.baseName.replaceAll(/\\.targetcoverage$/, '')
    """
    # 性别推断
    SEX_RESULT=\$(cnvkit.py sex ${target_coverage} ${antitarget_coverage} 2>&1)

    # 解析性别结果
    INFERRED_SEX=\$(echo "\${SEX_RESULT}" | grep -oP '(Male|Female)' | head -1 || echo "Unknown")

    # 计算 X/Y 染色体覆盖度比值
    X_COV=\$(grep -E "^(chr)?X\\s" ${target_coverage} 2>/dev/null | awk '{sum+=\$5; count++} END {if(count>0) print sum/count; else print 0}')
    Y_COV=\$(grep -E "^(chr)?Y\\s" ${target_coverage} 2>/dev/null | awk '{sum+=\$5; count++} END {if(count>0) print sum/count; else print 0}')

    # 输出 JSON 格式结果
    cat <<EOF > ${sample_id}.sex.json
    {
        "sample_id": "${sample_id}",
        "inferred_sex": "\${INFERRED_SEX}",
        "x_coverage": \${X_COV:-0},
        "y_coverage": \${Y_COV:-0},
        "raw_output": "\${SEX_RESULT}"
    }
    EOF
    """
}

