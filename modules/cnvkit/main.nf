// CNVKIT 模块
// 用途：基于测序深度的拷贝数变异 (CNV) 检测
// 工具：CNVkit
//
// 分析流程 (建立基线)：
//   1. CNVKIT_REFERENCE_TARGET - 从用户 BED 生成 target BED
//   2. CNVKIT_REFERENCE_ANTITARGET - 从 target BED 生成 antitarget BED
//   3. CNVKIT_REFERENCE_COVERAGE - 批量计算覆盖度 (每个样本并行)
//   4. CNVKIT_REFERENCE - 合并覆盖度生成基线文件
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
//   - 建立基线需要一批正常样本 (建议 5-10 个以上)
//   - 用户 BED 文件会自动转换为 target/antitarget 格式
//   - 支持自动性别检测和性别特异分析

// ============================================================================
// CNVKIT_REFERENCE - 建立基线
// ============================================================================
process CNVKIT_REFERENCE {
    tag "CNVKIT_REFERENCE"
    label 'process_high'
    label 'cnvkit'

    input:
        path alignments          // 批次 BAM/CRAM 比对文件
        path alignment_indices   // 比对文件索引
        path user_bed            // 用户提供的捕获区域 BED 文件
        path access_bed          // 可访问区域 BED (可选)
        path fasta               // 参考基因组
        val annotate             // 是否添加基因注释 (boolean)
        val split_size           // 分割大小 (默认 5000)
        val min_target_size      // 最小目标大小 (默认 10000)
        val is_female_reference  // 是否使用女性参考 (默认 true)
        val output_dir           // 输出目录

    output:
        path "reference.cnn", emit: reference

    when:
    output_dir != 'NO_OUTPUT'

    script:
    def min_size = min_target_size ?: 10000
    def split = split_size ?: 5000
    def annotate_cmd = annotate ? '--annotate' : ''
    def access_cmd = access_bed.name != 'NO_FILE' ? "--access ${access_bed}" : ''
    def female_cmd = is_female_reference ? '--female-reference' : ''
    """
    # 生成 target BED
    cnvkit.py target \\
        ${user_bed} \\
        --output targets.bed \\
        --split ${split} \\
        ${annotate_cmd}

    # 生成 antitarget BED
    cnvkit.py antitarget \\
        targets.bed \\
        --output antitargets.bed \\
        --min-target-size ${min_size} \\
        ${access_cmd}

    # 批量计算覆盖度
    mkdir -p coverage_files

    for aln in ${alignments}; do
        sample_id=\$(basename "\${aln}" | sed 's/\\.(marked|bam|cram)$//')

        if [[ "\${aln}" == *.cram ]]; then
            ref_cmd="--fasta ${fasta}"
        else
            ref_cmd=""
        fi

        cnvkit.py coverage \\
            "\${aln}" \\
            targets.bed \\
            \${ref_cmd} \\
            --processes ${task.cpus} \\
            --output coverage_files/\${sample_id}.targetcoverage.cnn

        cnvkit.py coverage \\
            "\${aln}" \\
            antitargets.bed \\
            \${ref_cmd} \\
            --processes ${task.cpus} \\
            --output coverage_files/\${sample_id}.antitargetcoverage.cnn
    done

    # 合并生成参考基线
    cnvkit.py reference \\
        coverage_files/*.targetcoverage.cnn \\
        coverage_files/*.antitargetcoverage.cnn \\
        --targets targets.bed \\
        --antitargets antitargets.bed \\
        ${female_cmd} \\
        --output reference.cnn
    """
}

// ============================================================================
// process CNVKIT_BATCH
// ============================================================================
process CNVKIT_BATCH {
    tag "CNVKIT_BATCH on ${alignment.baseName}"
    label 'process_low'
    label 'cnvkit'

    input:
        path alignment          // BAM/CRAM 比对文件
        path alignment_index    // 比对文件索引
        path user_bed           // 用户提供的捕获区域 BED 文件
        val annotate            // 是否添加基因注释 (boolean)
        val split_size          // 分割大小 (默认 5000)
        path access_bed         // 可访问区域 BED (可选，用于限制 antitarget 范围)
        val min_target_size     // 最小目标大小 (默认 10000)
        path fasta              // 参考基因组 (CRAM 输入时需要)
        path reference          // 基线文件 (reference.cnn)
        val method              // 分段方法: 'cbs' (默认), 'haar', 'flasso'
        val threshold           // 分段阈值 (默认自动)
        val ploidy              // 倍性 (默认 2)
        val is_male             // 样本是否为男性 (boolean)
        val drop_outliers       // 是否丢弃离群值 (boolean, 默认 true)
        val output_dir          // 输出目录

    output:
        path "targets.bed", emit: target_bed
        path "antitargets.bed", emit: antitarget_bed
        path "${sample_id}.targetcoverage.cnn", emit: target_coverage
        path "${sample_id}.antitargetcoverage.cnn", emit: antitarget_coverage
        path "${sample_id}.cnr", emit: cnr
        path "${sample_id}.exon.call.cns", emit: exon_call_cns
        path "${sample_id}.seg.cns", emit: seg_cns
        path "${sample_id}.seg.call.cns", emit: seg_call_cns

    when:
    output_dir != 'NO_OUTPUT'

    script:
    def min_size = min_target_size ?: 10000
    def access_cmd = access_bed.name != 'NO_FILE' ? "--access ${access_bed}" : ''
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam|cram)$/, '')
    def ref_cmd = alignment.name.endsWith('.cram') ? "--fasta ${fasta}" : ''
    def seg_method = method ?: 'cbs'
    def thresh_cmd = threshold ? "--threshold ${threshold}" : ''
    def male_cmd = is_male ? '--male-reference' : ''
    def drop_cmd = drop_outliers ? '--drop-low-coverage' : ''
    def annotate_cmd = annotate ? '--annotate' : ''
    def split = split_size ?: 5000
    """
    cnvkit.py target \\
        ${user_bed} \\
        --output targets.bed \\
        --split ${split} \\
        ${annotate_cmd}

    cnvkit.py antitarget \\
        targets.bed \\
        --output antitargets.bed \\
        --min-target-size ${min_size} \\
        ${access_cmd}

    cnvkit.py coverage \\
        ${alignment} \\
        targets.bed \\
        ${ref_cmd} \\
        --processes ${task.cpus} \\
        --output ${sample_id}.targetcoverage.cnn

    cnvkit.py coverage \\
        ${alignment} \\
        antitargets.bed \\
        ${ref_cmd} \\
        --processes ${task.cpus} \\
        --output ${sample_id}.antitargetcoverage.cnn

    cnvkit.py fix \\
        ${sample_id}.targetcoverage.cnn \\
        ${sample_id}.antitargetcoverage.cnn \\
        ${reference} \\
        --output ${sample_id}.cnr

    cnvkit.py call \\
        ${sample_id}.cnr \\
        --ploidy 2 \\
        ${male_cmd} \\
        ${drop_cmd} \\
        --output ${sample_id}.exon.call.cns

    cnvkit.py segment \\
        ${sample_id}.cnr \\
        --method ${seg_method} \\
        --processes ${task.cpus} \\
        ${thresh_cmd} \\
        --output ${sample_id}.seg.cns

    cnvkit.py call \\
        ${sample_id}.seg.cns \\
        --ploidy 2 \\
        ${male_cmd} \\
        ${drop_cmd} \\
        --output ${sample_id}.seg.call.cns
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
        val output_dir          // 输出目录

    output:
        path "*.sex.json", emit: sex_json

    when:
    output_dir != 'NO_OUTPUT'

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

