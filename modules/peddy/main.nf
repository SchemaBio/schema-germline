// PEDDY 模块
// 用途：亲缘关系验证和样本质量控制
// 工具：Peddy
// 输入：gVCF/VCF 文件 + PED 家系文件 + 参考基因组
// 输出：亲缘关系验证报告 + 性别验证 + 样本 QC 报告
//
// 说明：
//   - Peddy 用于验证样本之间的亲缘关系（亲子、兄弟姐妹等）
//   - 检测样本性别是否与 PED 文件一致
//   - 检测样本异常（样本交换、污染、意外亲缘关系等）
//   - 基于 VCF 中变异数据计算样本间的 IBS (Identity by State)
//   - 支持批量处理多个样本的家系分析

// ============================================================================
// PEDDY_CHECK - 亲缘关系验证
// ============================================================================
process PEDDY_CHECK {
    tag "PEDDY on ${ped.baseName}"
    label 'process_medium'
    label 'peddy'

    input:
        path vcfs                 // 多样本 VCF/gVCF 文件列表
        path vcfs_tbi             // VCF 索引文件列表
        path ped                  // PED 家系文件 (PLINK 格式)
        path fasta                // 参考基因组 FASTA
        path fasta_fai            // 参考基因组索引 (.fai)
        val genome_assembly        // 基因组版本: GRCh37 或 GRCh38
        val output_dir             // 输出目录 (可选)

    output:
        path "*.peddy.*.csv", emit: csv_results
        path "*.peddy.html", emit: html_report
        path "*.peddy-summary.json", emit: summary_json
        path "*.peddy.ped", emit: updated_ped

    script:
    def assembly = genome_assembly ?: 'GRCh38'
    """
    # 创建临时目录存放 VCF 文件
    mkdir -p vcf_temp

    # 复制 VCF 文件到临时目录
    # Peddy 需要所有样本的 VCF 文件在同一个目录
    for vcf in ${vcfs}; do
        cp \$vcf vcf_temp/
    done
    for tbi in ${vcfs_tbi}; do
        cp \$tbi vcf_temp/
    done

    # 运行 Peddy 亲缘关系验证
    # --sites: 使用 Peddy 内置的 SNP 位点或自定义位点文件
    # --plot: 生成可视化 HTML 报告
    # Peddy 会自动检测:
    #   - 亲缘关系是否与 PED 文件一致
    #   - 性别是否与 PED 文件一致
    #   - 样本间 IBS 矩阵
    #   - 异常亲缘关系或样本交换
    peddy \\
        --vcfs vcf_temp/*.vcf.gz \\
        --ped ${ped} \\
        --fasta ${fasta} \\
        --plot \\
        --out-dir . \\
        --cores ${task.cpus}

    # 清理临时目录
    rm -rf vcf_temp

    # 输出文件说明:
    #   *.peddy.ped.csv: 更新后的 PED 文件，包含验证结果
    #   *.peddy.relations.csv: 亲缘关系验证结果
    #   *.peddy.sex.csv: 性别验证结果
    #   *.peddy.samples.csv: 样本 QC 信息
    #   *.peddy.html: 可视化 HTML 报告
    #   *.peddy-summary.json: JSON 格式摘要
    #   *.peddy.ped: 更新后的 PED 文件 (如有修正)
    """
}


// ============================================================================
// PEDDY_PAIRWISE - 两样本亲缘关系分析
// ============================================================================
process PEDDY_PAIRWISE {
    tag "PEDDY_PAIRWISE: ${sample1} vs ${sample2}"
    label 'process_low'
    label 'peddy'

    input:
        val sample1               // 样本1 ID
        path vcf1                 // 样本1 VCF 文件
        path vcf1_tbi             // 样本1 VCF 索引
        val sample2               // 样本2 ID
        path vcf2                 // 样本2 VCF 文件
        path vcf2_tbi             // 样本2 VCF 索引
        val expected_relation      // 预期亲缘关系: parent-child, siblings, unrelated, etc.
        path fasta                // 参考基因组 FASTA
        path fasta_fai            // 参考基因组索引 (.fai)
        val output_dir             // 输出目录 (可选)

    output:
        path "*pairwise.csv", emit: pairwise_result
        path "*pairwise.json", emit: pairwise_json

    script:
    """
    # 创建临时目录
    mkdir -p vcf_temp

    # 复制 VCF 文件并重命名
    cp ${vcf1} vcf_temp/${sample1}.vcf.gz
    cp ${vcf1_tbi} vcf_temp/${sample1}.vcf.gz.tbi
    cp ${vcf2} vcf_temp/${sample2}.vcf.gz
    cp ${vcf2_tbi} vcf_temp/${sample2}.vcf.gz.tbi

    # 创建临时 PED 文件
    # PED 格式: FamilyID SampleID FatherID MotherID Sex Phenotype
    # 对于两样本分析，假设为同一家系
    cat > temp.ped << EOF
    F1 ${sample1} 0 0 0 0
    F1 ${sample2} 0 0 0 0
    EOF

    # 运行 Peddy
    peddy \\
        --vcfs vcf_temp/*.vcf.gz \\
        --ped temp.ped \\
        --fasta ${fasta} \\
        --cores ${task.cpus}

    # 提取两样本间的亲缘关系结果
    # 计算实际亲缘关系和预期关系是否一致
    python3 << 'PYEOF'
    import json
    import csv
    import os

    # 读取 peddy 输出的关系文件
    relations_file = "F1.peddy.relations.csv"
    if os.path.exists(relations_file):
        with open(relations_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if (row['sample_a'] == '${sample1}' and row['sample_b'] == '${sample2}') or \
                   (row['sample_a'] == '${sample2}' and row['sample_b'] == '${sample1}'):
                    result = {
                        'sample1': '${sample1}',
                        'sample2': '${sample2}',
                        'ibs0': float(row.get('ibs0', 0)),
                        'ibs1': float(row.get('ibs1', 0)),
                        'ibs2': float(row.get('ibs2', 0)),
                        'rel_type': row.get('rel_type', 'unknown'),
                        'expected_relation': '${expected_relation}',
                        'verified': row.get('rel_type', '').replace('_', '-') == '${expected_relation}' or \
                                    '${expected_relation}' == 'any'
                    }

                    # 写入 CSV
                    with open('${sample1}_${sample2}_pairwise.csv', 'w', newline='') as out:
                        writer = csv.DictWriter(out, fieldnames=result.keys())
                        writer.writeheader()
                        writer.writerow(result)

                    # 写入 JSON
                    with open('${sample1}_${sample2}_pairwise.json', 'w') as out:
                        json.dump(result, out, indent=2)
                    break
    PYEOF

    # 清理
    rm -rf vcf_temp temp.ped F1.peddy.*
    """
}


// ============================================================================
// PEDDY_FAMILY - 家系亲缘关系分析
// ============================================================================
process PEDDY_FAMILY {
    tag "PEDDY_FAMILY: ${family_id}"
    label 'process_medium'
    label 'peddy'

    input:
        val family_id              // 家系 ID
        path vcfs                  // 家系所有样本的 VCF 文件列表
        path vcfs_tbi              // VCF 索引文件列表
        val sample_ids             // 样本 ID 列表 (逗号分隔)
        val father_id              // 父样本 ID
        val mother_id              // 母样本 ID
        val proband_ids            // 先证者/子代样本 ID 列表 (逗号分隔)
        path fasta                 // 参考基因组 FASTA
        path fasta_fai             // 参考基因组索引 (.fai)
        val output_dir              // 输出目录 (可选)

    output:
        path "*family.csv", emit: family_result
        path "*family.json", emit: family_json
        path "*family.html", emit: family_report

    script:
    """
    # 创建临时目录
    mkdir -p vcf_temp

    # 复制 VCF 文件到临时目录
    for vcf in ${vcfs}; do
        sample_name=\$(basename \$vcf .vcf.gz)
        cp \$vcf vcf_temp/\$sample_name.vcf.gz
    done
    for tbi in ${vcfs_tbi}; do
        cp \$tbi vcf_temp/
    done

    # 创建家系 PED 文件
    # PED 格式: FamilyID SampleID FatherID MotherID Sex Phenotype
    # Sex: 1=male, 2=female, 0=unknown
    # Phenotype: 1=unaffected, 2=affected, 0=unknown
    cat > ${family_id}.ped << EOF
    ${family_id} ${father_id} 0 0 1 0
    ${family_id} ${mother_id} 0 0 2 0
    EOF

    # 添加先证者/子代样本
    # 默认性别为未知 (0)，可根据实际调整
    for proband in ${proband_ids}; do
        echo "${family_id} \$proband ${father_id} ${mother_id} 0 2" >> ${family_id}.ped
    done

    # 运行 Peddy
    peddy \\
        --vcfs vcf_temp/*.vcf.gz \\
        --ped ${family_id}.ped \\
        --fasta ${fasta} \\
        --plot \\
        --cores ${task.cpus}

    # 提取家系结果
    mv ${family_id}.peddy.html ${family_id}_family.html || true
    mv ${family_id}.peddy.relations.csv ${family_id}_family.csv || true
    mv ${family_id}.peddy-summary.json ${family_id}_family.json || true

    # 清理临时目录
    rm -rf vcf_temp
    """
}


// ============================================================================
// PEDDY_SUMMARY - 亲缘关系结果汇总
// ============================================================================
process PEDDY_SUMMARY {
    tag "PEDDY_SUMMARY"
    label 'process_low'
    label 'peddy'

    input:
        path csv_results           // Peddy 输出的 CSV 结果文件
        path json_results          // Peddy 输出的 JSON 结果文件

    output:
        path "peddy_final_report.csv", emit: final_csv
        path "peddy_final_report.json", emit: final_json
        path "peddy_issues.txt", emit: issues_report

    script:
    """
    # 汇总所有 Peddy 结果
    # 检查是否有异常样本或验证失败的亲缘关系

    echo "PEDDY Summary Report" > peddy_final_report.csv
    echo "Generated: \$(date)" >> peddy_final_report.csv
    echo "" >> peddy_final_report.csv

    # 合并所有 CSV 结果
    for csv in ${csv_results}; do
        if [ -f "\$csv" ]; then
            echo "=== \$csv ===" >> peddy_final_report.csv
            cat \$csv >> peddy_final_report.csv
            echo "" >> peddy_final_report.csv
        fi
    done

    # 合并所有 JSON 结果
    python3 << 'PYEOF'
    import json
    import glob
    import os

    all_results = []
    issues = []

    for json_file in glob.glob("*.peddy-summary.json"):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
                all_results.append(data)

                # 检查是否有问题
                if 'problems' in data:
                    for problem in data['problems']:
                        issues.append(f"{json_file}: {problem}")
        except Exception as e:
            issues.append(f"{json_file}: Error reading file - {str(e)}")

    # 写入合并的 JSON
    with open('peddy_final_report.json', 'w') as out:
        json.dump({
            'results': all_results,
            'total_samples': len(all_results),
            'issues_count': len(issues)
        }, out, indent=2)

    # 写入问题报告
    with open('peddy_issues.txt', 'w') as out:
        if issues:
            out.write("PEDDY Issues Report\\n")
            out.write("=" * 50 + "\\n\\n")
            for issue in issues:
                out.write(f"{issue}\\n")
        else:
            out.write("No issues detected.\\n")
            out.write("All samples passed pedigree verification.\\n")
    PYEOF
    """
}