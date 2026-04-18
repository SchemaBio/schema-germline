// CNVKIT 模块
// 用途：基于测序深度的拷贝数变异 (CNV) 检测
// 工具：CNVkit
//
// 分析流程 (建立基线)：
//   1. CNVKIT_TARGET_ANTITARGET - 从用户 BED 生成 target/antitarget BED
//   2. CNVKIT_COVERAGE - 计算覆盖度 (每个样本并行)
//   3. CNVKIT_REFERENCE_BUILD - 合并覆盖度生成基线文件
//
// 分析流程 (单样本模式，已有基线文件)：
//   1. CNVKIT_COVERAGE - 计算覆盖度
//   2. CNVKIT_FIX - 修正覆盖度偏差 (需要基线文件)
//   3. CNVKIT_SEGMENT - CBS 分段
//   4. CNVKIT_CALL - CNV 调用
//   5. CNVKIT_EXPORT_VCF - 导出 VCF 格式
//
// 说明：
//   - 建立基线需要一批正常样本 (建议 5-10 个以上)
//   - 用户 BED 文件会自动转换为 target/antitarget 格式
//   - 支持自动性别检测和性别特异分析

// ============================================================================
// CNVKIT_TARGET_ANTITARGET - 生成 target 和 antitarget BED 文件
// ============================================================================
process CNVKIT_TARGET_ANTITARGET {
    tag "CNVKIT_TARGET_ANTITARGET"
    label 'process_low'
    label 'cnvkit'

    input:
        path user_bed           // 用户提供的捕获区域 BED 文件
        val has_access_bed      // 是否有 access_bed 文件
        path access_bed         // 可访问区域 BED (可选)
        val annotate            // 是否添加基因注释 (boolean)
        val split_size          // 分割大小 (默认 5000)
        val min_target_size     // 最小目标大小 (默认 10000)

    output:
        path "targets.bed", emit: target_bed
        path "antitargets.bed", emit: antitarget_bed

    script:
    def min_size = min_target_size ?: 10000
    def split = split_size ?: 5000
    def annotate_cmd = annotate ? '--annotate' : ''
    """
    # 生成 target BED
    cnvkit.py target \\
        ${user_bed} \\
        --output targets.bed \\
        --split ${split} \\
        ${annotate_cmd}

    # 生成 antitarget BED
    if [ "${has_access_bed}" == "true" ] && [ -f "${access_bed}" ]; then
        cnvkit.py antitarget \\
            targets.bed \\
            --output antitargets.bed \\
            --min-target-size ${min_size} \\
            --access ${access_bed}
    else
        cnvkit.py antitarget \\
            targets.bed \\
            --output antitargets.bed \\
            --min-target-size ${min_size}
    fi
    """
}

// ============================================================================
// CNVKIT_COVERAGE - 计算单个样本的覆盖度
// ============================================================================
process CNVKIT_COVERAGE {
    tag "CNVKIT_COVERAGE on ${alignment.baseName}"
    label 'process_medium'
    label 'cnvkit'

    input:
        path alignment          // BAM 比对文件
        path alignment_index    // 比对文件索引
        path target_bed         // target BED 文件
        path antitarget_bed     // antitarget BED 文件

    output:
        path "${sample_id}.targetcoverage.cnn", emit: target_coverage
        path "${sample_id}.antitargetcoverage.cnn", emit: antitarget_coverage

    script:
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam)$/, '')
    """
    # 计算目标区域覆盖度
    cnvkit.py coverage \\
        ${alignment} \\
        ${target_bed} \\
        --processes ${task.cpus} \\
        --output ${sample_id}.targetcoverage.cnn

    # 计算反目标区域覆盖度
    cnvkit.py coverage \\
        ${alignment} \\
        ${antitarget_bed} \\
        --processes ${task.cpus} \\
        --output ${sample_id}.antitargetcoverage.cnn
    """
}

// ============================================================================
// CNVKIT_REFERENCE_BUILD - 合并覆盖度文件生成基线
// ============================================================================
process CNVKIT_REFERENCE_BUILD {
    tag "CNVKIT_REFERENCE_BUILD: ${target_coverages.size()} samples"
    label 'process_high'
    label 'cnvkit'

    input:
        path target_coverages   // 多个 target coverage 文件
        path antitarget_coverages // 多个 antitarget coverage 文件
        path target_bed         // target BED 文件
        path antitarget_bed     // antitarget BED 文件
        val is_female_reference // 是否使用女性参考 (默认 true)

    output:
        path "reference.cnn", emit: reference

    script:
    def female_cmd = is_female_reference ? '--female-reference' : ''
    """
    # 合并生成参考基线
    cnvkit.py reference \\
        ${target_coverages} \\
        ${antitarget_coverages} \\
        --targets ${target_bed} \\
        --antitargets ${antitarget_bed} \\
        ${female_cmd} \\
        --output reference.cnn
    """
}

// ============================================================================
// CNVKIT_FIX - 修正覆盖度偏差
// ============================================================================
process CNVKIT_FIX {
    tag "CNVKIT_FIX on ${alignment.baseName}"
    label 'process_low'
    label 'cnvkit'

    input:
        path alignment          // BAM 比对文件 (用于获取 sample_id)
        path target_coverage    // 目标区域覆盖度文件
        path antitarget_coverage // 反目标区域覆盖度文件
        path reference          // 基线文件 (reference.cnn)

    output:
        path "${sample_id}.cnr", emit: cnr

    script:
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam)$/, '')
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

    input:
        path cnr                // 修正后的覆盖度比率文件
        val method              // 分段方法: 'cbs' (默认), 'haar', 'flasso'
        val threshold           // 分段阈值 (可选)

    output:
        path "${sample_id}.seg.cns", emit: seg_cns

    script:
    def sample_id = cnr.baseName.replaceAll(/\.cnr$/, '')
    def seg_method = method ?: 'cbs'
    def thresh_cmd = threshold ? "--threshold ${threshold}" : ''
    """
    cnvkit.py segment \\
        ${cnr} \\
        --method ${seg_method} \\
        --processes ${task.cpus} \\
        ${thresh_cmd} \\
        --output ${sample_id}.seg.cns
    """
}

// ============================================================================
// CNVKIT_CALL - CNV 调用
// ============================================================================
process CNVKIT_CALL {
    tag "CNVKIT_CALL on ${cns.baseName}"
    label 'process_low'
    label 'cnvkit'

    input:
        path cns                // 分段文件 (.seg.cns 或 .cns)
        val ploidy              // 倍性 (默认 2)
        val is_male             // 样本是否为男性 (boolean)
        val drop_outliers       // 是否丢弃离群值 (boolean, 默认 true)

    output:
        path "${sample_id}.call.cns", emit: call_cns

    script:
    def sample_id = cns.baseName.replaceAll(/\.seg\.cns$/, '').replaceAll(/\.cns$/, '')
    def male_cmd = is_male ? '--male-reference' : ''
    def drop_cmd = drop_outliers ? '--drop-low-coverage' : ''
    def ploid = ploidy ?: 2
    """
    cnvkit.py call \\
        ${cns} \\
        --ploidy ${ploid} \\
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

    input:
        path call_cns           // CNV 调用结果
        path cnr                // 修正后的覆盖度比率文件

    output:
        path "${sample_id}.cnv.vcf.gz", emit: vcf
        path "${sample_id}.cnv.vcf.gz.tbi", emit: vcf_tbi

    script:
    def sample_id = call_cns.baseName.replaceAll(/\.call\.cns$/, '')
    """
    cnvkit.py export vcf \\
        ${call_cns} \\
        --sample-id ${sample_id} \\
        -o ${sample_id}.cnv.vcf

    bgzip -c ${sample_id}.cnv.vcf > ${sample_id}.cnv.vcf.gz
    tabix -p vcf ${sample_id}.cnv.vcf.gz
    """
}

// ============================================================================
// CNVKIT_SEX - 性别检测
// ============================================================================
process CNVKIT_SEX {
    tag "CNVKIT_SEX on ${target_coverage.baseName}"
    label 'process_low'
    label 'cnvkit'

    input:
        path target_coverage    // 目标区域覆盖度文件
        path antitarget_coverage // 反目标区域覆盖度文件

    output:
        path "*.sex.json", emit: sex_json

    script:
    def sample_id = target_coverage.baseName.replaceAll(/\.targetcoverage$/, '')
    """
    # 性别推断
    SEX_RESULT=$(cnvkit.py sex ${target_coverage} ${antitarget_coverage} 2>&1)

    # 解析性别结果
    INFERRED_SEX=$(echo "${SEX_RESULT}" | grep -oP '(Male|Female)' | head -1 || echo "Unknown")

    # 计算 X/Y 染色体覆盖度比值
    X_COV=$(grep -E "^(chr)?X\\s" ${target_coverage} 2>/dev/null | awk '{sum+=$5; count++} END {if(count>0) print sum/count; else print 0}')
    Y_COV=$(grep -E "^(chr)?Y\\s" ${target_coverage} 2>/dev/null | awk '{sum+=$5; count++} END {if(count>0) print sum/count; else print 0}')

    # 输出 JSON 格式结果
    cat <<EOF > ${sample_id}.sex.json
    {
        "sample_id": "${sample_id}",
        "inferred_sex": "${INFERRED_SEX}",
        "x_coverage": ${X_COV:-0},
        "y_coverage": ${Y_COV:-0},
        "raw_output": "${SEX_RESULT}"
    }
    EOF
    """
}

// ============================================================================
// CNVKIT_BATCH - 单样本完整流程 (已有基线)
// ============================================================================
process CNVKIT_BATCH {
    tag "CNVKIT_BATCH on ${alignment.baseName}"
    label 'process_medium'
    label 'cnvkit'

    input:
        path alignment          // BAM 比对文件
        path alignment_index    // 比对文件索引
        path target_bed         // target BED 文件
        path antitarget_bed     // antitarget BED 文件
        path reference          // 基线文件 (reference.cnn)
        val method              // 分段方法: 'cbs' (默认), 'haar', 'flasso'
        val threshold           // 分段阈值 (可选)
        val ploidy              // 倍性 (默认 2)
        val is_male             // 样本是否为男性 (boolean)
        val drop_outliers       // 是否丢弃离群值 (boolean, 默认 true)
        val output_dir          // 输出目录

    output:
        path "${sample_id}.targetcoverage.cnn", emit: target_coverage
        path "${sample_id}.antitargetcoverage.cnn", emit: antitarget_coverage
        path "${sample_id}.cnr", emit: cnr
        path "${sample_id}.seg.cns", emit: seg_cns
        path "${sample_id}.seg.call.cns", emit: call_cns
        path "${sample_id}.cnv.vcf.gz", emit: vcf
        path "${sample_id}.cnv.vcf.gz.tbi", emit: vcf_tbi

    when:
    output_dir != 'NO_OUTPUT'

    script:
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam)$/, '')
    def seg_method = method ?: 'cbs'
    def thresh_cmd = threshold ? "--threshold ${threshold}" : ''
    def male_cmd = is_male ? '--male-reference' : ''
    def drop_cmd = drop_outliers ? '--drop-low-coverage' : ''
    def ploid = ploidy ?: 2
    """
    # 计算覆盖度
    cnvkit.py coverage \\
        ${alignment} \\
        ${target_bed} \\
        --processes ${task.cpus} \\
        --output ${sample_id}.targetcoverage.cnn

    cnvkit.py coverage \\
        ${alignment} \\
        ${antitarget_bed} \\
        --processes ${task.cpus} \\
        --output ${sample_id}.antitargetcoverage.cnn

    # 修正覆盖度偏差
    cnvkit.py fix \\
        ${sample_id}.targetcoverage.cnn \\
        ${sample_id}.antitargetcoverage.cnn \\
        ${reference} \\
        --output ${sample_id}.cnr

    # CBS 分段
    cnvkit.py segment \\
        ${sample_id}.cnr \\
        --method ${seg_method} \\
        --processes ${task.cpus} \\
        ${thresh_cmd} \\
        --output ${sample_id}.seg.cns

    # CNV 调用
    cnvkit.py call \\
        ${sample_id}.seg.cns \\
        --ploidy ${ploid} \\
        ${male_cmd} \\
        ${drop_cmd} \\
        --output ${sample_id}.seg.call.cns

    # 导出 VCF
    cnvkit.py export vcf \\
        ${sample_id}.seg.call.cns \\
        --sample-id ${sample_id} \\
        -o ${sample_id}.cnv.vcf

    bgzip -c ${sample_id}.cnv.vcf > ${sample_id}.cnv.vcf.gz
    tabix -p vcf ${sample_id}.cnv.vcf.gz
    """
}