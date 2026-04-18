// GLnexus 模块
// 用途：合并多个 gVCF 文件进行 cohort joint calling
// 包含的 process：
//   GLNEXUS_JOINTCALL - 合并多个 g.vcf.gz 文件生成 cohort VCF

process GLNEXUS_JOINTCALL {
    tag "GLNEXUS_JOINTCALL: ${gvcfs.size()} samples"
    label 'glnexus'
    label 'process_high'

    input:
        path gvcfs           // 多个 g.vcf.gz 文件 (集合或列表)
        path gvcfs_indices   // 对应的 .tbi 索引文件
        val config_mode      // GLnexus 配置模式: 'deepvariant' 或 'gatk' 等
        val output_prefix    // 输出文件前缀

    output:
        path "${output_prefix}.vcf.gz", emit: vcf
        path "${output_prefix}.vcf.gz.tbi", emit: vcf_tbi

    script:
    """
    # 创建临时目录存放 gVCF 文件
    mkdir -p gvcf_input

    # 复制所有 gVCF 文件到临时目录
    # GLnexus 需要所有文件在同一目录下以便正确处理
    for gvcf in ${gvcfs}; do
        cp \$gvcf gvcf_input/
    done

    # 复制索引文件
    for tbi in ${gvcfs_indices}; do
        cp \$tbi gvcf_input/
    done

    # 运行 GLnexus joint calling
    # --config: 指定配置模式 (deepvariant/gatk/gatk4_gvcf 等)
    # --dir: 从目录读取所有 gVCF 文件
    # --ped: 可选的 PED 文件用于 trio 分析
    glnexus_cli \\
        --config ${config_mode} \\
        --dir gvcf_input \\
        --threads ${task.cpus} \\
        ${output_prefix}.vcf.gz

    # 创建 VCF 索引
    tabix -p vcf ${output_prefix}.vcf.gz

    # 清理临时目录
    rm -rf gvcf_input
    """
}

// GLNEXUS_JOINTCALL_TRIO - Trio 家系联合变异检测
// 专门用于 trio 分析，可以推断 de novo 变异
process GLNEXUS_JOINTCALL_TRIO {
    tag "GLNEXUS_JOINTCALL_TRIO: ${ped_file.baseName}"
    label 'glnexus'
    label 'process_high'

    input:
        path gvcfs           // 多个 g.vcf.gz 文件 (proband + parents)
        path gvcfs_indices   // 对应的 .tbi 索引文件
        path ped_file        // PED 文件描述家系关系
        val config_mode      // GLnexus 配置模式
        val output_prefix    // 输出文件前缀

    output:
        path "${output_prefix}.vcf.gz", emit: vcf
        path "${output_prefix}.vcf.gz.tbi", emit: vcf_tbi
        path "${output_prefix}.de_novo.vcf.gz", emit: de_novo_vcf, optional: true
        path "${output_prefix}.de_novo.vcf.gz.tbi", emit: de_novo_tbi, optional: true

    script:
    """
    # 创建临时目录存放 gVCF 文件
    mkdir -p gvcf_input

    # 复制所有 gVCF 文件到临时目录
    for gvcf in ${gvcfs}; do
        cp \$gvcf gvcf_input/
    done

    # 复制索引文件
    for tbi in ${gvcfs_indices}; do
        cp \$tbi gvcf_input/
    done

    # 运行 GLnexus joint calling (trio 模式)
    # --ped: 指定 PED 文件用于 trio 分析
    # --multiparent: 支持多个家系同时分析
    glnexus_cli \\
        --config ${config_mode} \\
        --dir gvcf_input \\
        --ped ${ped_file} \\
        --threads ${task.cpus} \\
        ${output_prefix}.vcf.gz

    # 创建 VCF 索引
    tabix -p vcf ${output_prefix}.vcf.gz

    # 提取 de novo 变异 (可选)
    # 如果 VCF 中包含 de novo 标记 (由 GLnexus trio 模式产生)
    if bcftools view -h ${output_prefix}.vcf.gz | grep -q "de novo"; then
        bcftools view \\
            -i 'INFO/de_novo=1' \\
            -Oz \\
            -o ${output_prefix}.de_novo.vcf.gz \\
            ${output_prefix}.vcf.gz

        tabix -p vcf ${output_prefix}.de_novo.vcf.gz
    fi

    # 清理临时目录
    rm -rf gvcf_input
    """
}

// GLNEXUS_MERGE_BATCH - 大规模队列的批量合并
// 当样本数量较多时 (>100)，建议分批合并后再合并
process GLNEXUS_MERGE_BATCH {
    tag "GLNEXUS_MERGE_BATCH: batch ${batch_id}"
    label 'glnexus'
    label 'process_high'

    input:
        path gvcfs           // 当前批次的 g.vcf.gz 文件
        path gvcfs_indices   // 对应的 .tbi 索引文件
        val config_mode      // GLnexus 配置模式
        val batch_id         // 批次 ID

    output:
        path "batch_${batch_id}.vcf.gz", emit: batch_vcf
        path "batch_${batch_id}.vcf.gz.tbi", emit: batch_tbi

    script:
    """
    mkdir -p gvcf_input

    for gvcf in ${gvcfs}; do
        cp \$gvcf gvcf_input/
    done

    for tbi in ${gvcfs_indices}; do
        cp \$tbi gvcf_input/
    done

    glnexus_cli \\
        --config ${config_mode} \\
        --dir gvcf_input \\
        --threads ${task.cpus} \\
        batch_${batch_id}.vcf.gz

    tabix -p vcf batch_${batch_id}.vcf.gz

    rm -rf gvcf_input
    """
}

// GLNEXUS_FINAL_MERGE - 合并多个批次的中间 VCF
// 用于大规模队列的两阶段合并策略
process GLNEXUS_FINAL_MERGE {
    tag "GLNEXUS_FINAL_MERGE: ${batch_vcfs.size()} batches"
    label 'glnexus'
    label 'process_huge'

    input:
        path batch_vcfs      // 多个批次的 VCF 文件
        path batch_vcfs_tbi  // 对应的索引文件
        val output_prefix    // 最终输出文件前缀

    output:
        path "${output_prefix}.cohort.vcf.gz", emit: final_vcf
        path "${output_prefix}.cohort.vcf.gz.tbi", emit: final_tbi

    script:
    """
    # 使用 bcftools merge 合并多个批次的 VCF
    # 注意: 所有批次的 VCF 必须有相同的样本集合和位点集合
    # 这是一个简化的合并方法，实际大规模队列可能需要更复杂的策略

    # 合并 VCF 文件列表
    vcf_list=""
    for vcf in ${batch_vcfs}; do
        vcf_list="\$vcf_list \$vcf"
    done

    bcftools merge \\
        --threads ${task.cpus} \\
        -Oz \\
        -o ${output_prefix}.cohort.vcf.gz \\
        \$vcf_list

    tabix -p vcf ${output_prefix}.cohort.vcf.gz
    """
}