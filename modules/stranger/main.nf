// STRANGER 模块
// 用途：STR 变异注释
// 工具：Stranger
// 输入：ExpansionHunter 输出的 VCF 文件
// 输出：注释后的 STR VCF，包含疾病信息、正常/致病范围等
//
// 说明：
//   - Stranger 专门用于注释 ExpansionHunter 生成的 STR VCF
//   - 添加疾病名称 (Disease)、正常范围 (NormalMax)、致病范围 (PathologicMin)
//   - 使用 variant catalog 文件提供 STR 位点定义和疾病关联信息
//   - 支持自定义或使用 Stranger 内置的 variant catalog

// ============================================================================
// STRANGER_ANNOTATE - STR 注释
// ============================================================================
process STRANGER_ANNOTATE {
    tag "STRANGER on ${vcf.baseName}"
    label 'process_low'
    label 'cnvkit'  // Stranger 包含在 cnvkit 容器中

    input:
        path vcf                   // ExpansionHunter 输出的 VCF 文件
        path vcf_tbi               // VCF 索引
        path variant_catalog       // STR 位点定义文件 (JSON 格式，与 ExpansionHunter 相同)
        val genome_assembly        // 基因组版本: GRCh37 或 GRCh38
        val output_dir             // 输出目录

    output:
        path "*.stranger.vcf.gz", emit: vcf
        path "*.stranger.vcf.gz.tbi", emit: vcf_tbi
        path "*.stranger.summary.txt", emit: summary

    script:
    def assembly = genome_assembly ?: 'GRCh38'
    def sample_id = vcf.baseName.replaceAll(/\.(vcf|vcf\.gz)$/, '')
    """
    # 运行 Stranger 注释
    # Stranger 使用与 ExpansionHunter 相同的 variant catalog
    # 注释字段包括：
    #   - STR_NORMAL_MAX: 正常范围最大重复数
    #   - STR_PATHOLOGIC_MIN: 致病范围最小重复数
    #   - Disease: 相关疾病名称
    #   - HGNC: 基因符号
    stranger \\
        --input ${vcf} \\
        --catalog ${variant_catalog} \\
        --genome ${assembly} \\
        --output ${sample_id}.stranger.vcf

    # 压缩 VCF 并创建索引
    bgzip -c ${sample_id}.stranger.vcf > ${sample_id}.stranger.vcf.gz
    tabix -p vcf ${sample_id}.stranger.vcf.gz

    # 清理原始 VCF
    rm -f ${sample_id}.stranger.vcf

    # 生成 STR 注释摘要
    echo "STRANGER Annotation Summary for ${sample_id}" > ${sample_id}.stranger.summary.txt
    echo "Input VCF: ${vcf}" >> ${sample_id}.stranger.summary.txt
    echo "Output VCF: ${sample_id}.stranger.vcf.gz" >> ${sample_id}.stranger.summary.txt
    echo "Genome assembly: ${assembly}" >> ${sample_id}.stranger.summary.txt
    echo "" >> ${sample_id}.stranger.summary.txt
    echo "Total STR variants annotated:" >> ${sample_id}.stranger.summary.txt
    zcat ${sample_id}.stranger.vcf.gz | grep -v "^#" | wc -l >> ${sample_id}.stranger.summary.txt
    echo "" >> ${sample_id}.stranger.summary.txt
    echo "STR loci with disease associations:" >> ${sample_id}.stranger.summary.txt
    zcat ${sample_id}.stranger.vcf.gz | grep -v "^#" | grep "Disease=" | cut -f8 | grep -o "Disease=[^;]*" | sort | uniq -c >> ${sample_id}.stranger.summary.txt
    echo "" >> ${sample_id}.stranger.summary.txt
    echo "Potential pathogenic STR expansions:" >> ${sample_id}.stranger.summary.txt
    zcat ${sample_id}.stranger.vcf.gz | grep -v "^#" | awk -F';' '{
        for(i=1;i<=NF;i++) {
            if($i ~ /STR_PATHOLOGIC_MIN=/) {
                path_min = substr($i, 19)
            }
            if($i ~ /REPCN=/) {
                repcn = substr($i, 7)
            }
        }
        if (repcn >= path_min && path_min > 0) print $0
    }' | head -20 >> ${sample_id}.stranger.summary.txt
    """
}


// ============================================================================
// STRANGER_FILTER - 筛选致病性 STR 扩展
// ============================================================================
process STRANGER_FILTER {
    tag "STRANGER_FILTER on ${vcf.baseName}"
    label 'process_low'
    label 'cnvkit'

    input:
        path vcf                   // Stranger 注释后的 VCF 文件
        path vcf_tbi               // VCF 索引
        val filter_mode            // 过滤模式: 'pathogenic', 'borderline', 'all_disease'
        val output_dir             // 输出目录

    output:
        path "*.stranger.filtered.vcf.gz", emit: filtered_vcf
        path "*.stranger.filtered.vcf.gz.tbi", emit: filtered_vcf_tbi
        path "*.stranger.pathogenic.txt", emit: report

    script:
    def sample_id = vcf.baseName.replaceAll(/\.(stranger\.vcf|stranger\.vcf\.gz)$/, '')
    def mode = filter_mode ?: 'pathogenic'
    """
    # 根据过滤模式筛选 STR 变异
    # pathogenic: 重复数 >= 致病阈值
    # borderline: 重复数介于正常上限和致病阈值之间
    # all_disease: 所有有疾病关联的 STR 位点

    if [ "${mode}" = "pathogenic" ]; then
        # 筛选致病性扩展
        zcat ${vcf} | grep "^#" > header.vcf
        zcat ${vcf} | grep -v "^#" | awk -F';' '{
            path_min = 0
            repcn = 0
            for(i=1;i<=NF;i++) {
                if($i ~ /STR_PATHOLOGIC_MIN=/) {
                    path_min = substr($i, 19)
                }
                if($i ~ /REPCN=/) {
                    repcn = substr($i, 7)
                }
            }
            if (repcn >= path_min && path_min > 0) print $0
        }' >> header.vcf
        cat header.vcf | bgzip -c > ${sample_id}.stranger.filtered.vcf.gz
        tabix -p vcf ${sample_id}.stranger.filtered.vcf.gz
        rm -f header.vcf

        echo "Pathogenic STR expansions for ${sample_id}" > ${sample_id}.stranger.pathogenic.txt
        zcat ${sample_id}.stranger.filtered.vcf.gz | grep -v "^#" >> ${sample_id}.stranger.pathogenic.txt

    elif [ "${mode}" = "borderline" ]; then
        # 筛选边界区域扩展
        zcat ${vcf} | grep "^#" > header.vcf
        zcat ${vcf} | grep -v "^#" | awk -F';' '{
            normal_max = 0
            path_min = 999999
            repcn = 0
            for(i=1;i<=NF;i++) {
                if($i ~ /STR_NORMAL_MAX=/) {
                    normal_max = substr($i, 16)
                }
                if($i ~ /STR_PATHOLOGIC_MIN=/) {
                    path_min = substr($i, 19)
                }
                if($i ~ /REPCN=/) {
                    repcn = substr($i, 7)
                }
            }
            if (repcn > normal_max && repcn < path_min && normal_max > 0) print $0
        }' >> header.vcf
        cat header.vcf | bgzip -c > ${sample_id}.stranger.filtered.vcf.gz
        tabix -p vcf ${sample_id}.stranger.filtered.vcf.gz
        rm -f header.vcf

        echo "Borderline STR expansions for ${sample_id}" > ${sample_id}.stranger.pathogenic.txt
        zcat ${sample_id}.stranger.filtered.vcf.gz | grep -v "^#" >> ${sample_id}.stranger.pathogenic.txt

    elif [ "${mode}" = "all_disease" ]; then
        # 筛选所有有疾病关联的位点
        zcat ${vcf} | grep "^#" > header.vcf
        zcat ${vcf} | grep -v "^#" | grep "Disease=" >> header.vcf
        cat header.vcf | bgzip -c > ${sample_id}.stranger.filtered.vcf.gz
        tabix -p vcf ${sample_id}.stranger.filtered.vcf.gz
        rm -f header.vcf

        echo "All disease-associated STR loci for ${sample_id}" > ${sample_id}.stranger.pathogenic.txt
        zcat ${sample_id}.stranger.filtered.vcf.gz | grep -v "^#" >> ${sample_id}.stranger.pathogenic.txt

    else
        # 默认输出全部
        cp ${vcf} ${sample_id}.stranger.filtered.vcf.gz
        cp ${vcf_tbi} ${sample_id}.stranger.filtered.vcf.gz.tbi
        echo "All STR variants for ${sample_id}" > ${sample_id}.stranger.pathogenic.txt
        zcat ${sample_id}.stranger.filtered.vcf.gz | grep -v "^#" >> ${sample_id}.stranger.pathogenic.txt
    fi
    """
}