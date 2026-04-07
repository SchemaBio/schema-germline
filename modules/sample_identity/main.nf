/*
 * 样本识别模块
 * 
 * 功能：基于 SNP 指纹位点进行样本身份验证和混淆检测
 * 工具：bcftools
 * 
 * 说明：使用预定义的 SNP 位点集提取基因型，生成样本指纹用于：
 *       1. 样本身份验证（同一样本不同批次比对）
 *       2. 样本混淆检测（不同样本间相似度）
 *       3. 亲缘关系初筛
 */

process SAMPLE_FINGERPRINT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 样本 VCF
    path  snp_panel                             // 指纹位点 VCF (如 1000G 常用位点)

    output:
    tuple val(meta), path("*.fingerprint.vcf.gz"), path("*.fingerprint.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.fingerprint.txt")                                     , emit: genotypes
    path "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # 提取指纹位点的基因型
    bcftools view \\
        --regions-file ${snp_panel} \\
        --output-type z \\
        --output ${prefix}.fingerprint.vcf.gz \\
        ${vcf}

    tabix -p vcf ${prefix}.fingerprint.vcf.gz

    # 生成简化的基因型文本 (位点 + GT)
    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n' \\
        ${prefix}.fingerprint.vcf.gz \\
        > ${prefix}.fingerprint.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

/*
 * SAMPLE_IDENTITY_CHECK - 样本间相似度比对
 * 
 * 比较两个样本的指纹，计算一致性
 */
process SAMPLE_IDENTITY_CHECK {
    tag "${meta1.id}_vs_${meta2.id}"
    label 'process_low'

    input:
    tuple val(meta1), path(fp1), path(fp1_idx), val(meta2), path(fp2), path(fp2_idx)

    output:
    tuple val(meta1), val(meta2), path("*.identity.txt"), emit: result
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta1.id}_vs_${meta2.id}"
    def concordance_threshold = params.identity_concordance_threshold ?: 0.95
    """
    # 合并两个样本的指纹 VCF
    bcftools merge \\
        --output-type z \\
        --output merged.vcf.gz \\
        ${fp1} ${fp2}

    tabix -p vcf merged.vcf.gz

    # 提取共同位点的基因型
    bcftools query \\
        -f '%CHROM\\t%POS\\t[%GT\\t]\\n' \\
        merged.vcf.gz \\
        > genotypes.txt

    # 计算一致性
    awk -F'\\t' '
    BEGIN { match=0; total=0; discord=0 }
    {
        gt1 = \$3
        gt2 = \$4
        # 跳过缺失基因型
        if (gt1 ~ /\\./ || gt2 ~ /\\./) next
        total++
        if (gt1 == gt2) {
            match++
        } else {
            discord++
        }
    }
    END {
        if (total > 0) {
            concordance = match / total
        } else {
            concordance = 0
        }
        printf "sample1: ${meta1.id}\\n"
        printf "sample2: ${meta2.id}\\n"
        printf "total_sites: %d\\n", total
        printf "concordant: %d\\n", match
        printf "discordant: %d\\n", discord
        printf "concordance: %.4f\\n", concordance
        if (concordance >= ${concordance_threshold}) {
            printf "status: SAME_SAMPLE\\n"
        } else if (concordance >= 0.4) {
            printf "status: RELATED\\n"
        } else {
            printf "status: DIFFERENT\\n"
        }
    }
    ' genotypes.txt > ${prefix}.identity.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

/*
 * SAMPLE_IDENTITY_MATRIX - 多样本相似度矩阵
 * 
 * 生成所有样本间的相似度矩阵，用于批量检测
 */
process SAMPLE_IDENTITY_MATRIX {
    tag "matrix"
    label 'process_medium'

    input:
    path fingerprint_vcfs  // 多个 fingerprint.vcf.gz 文件
    path fingerprint_tbis  // 对应的索引文件

    output:
    path "identity_matrix.tsv"  , emit: matrix
    path "identity_warnings.txt", emit: warnings
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def concordance_threshold = params.identity_concordance_threshold ?: 0.95
    """
    # 合并所有样本
    bcftools merge \\
        --output-type z \\
        --output all_samples.vcf.gz \\
        ${fingerprint_vcfs}

    tabix -p vcf all_samples.vcf.gz

    # 获取样本列表
    bcftools query -l all_samples.vcf.gz > samples.txt

    # 使用 bcftools gtcheck 计算样本间差异
    bcftools gtcheck \\
        --genotypes all_samples.vcf.gz \\
        --output gtcheck_output.txt

    # 解析 gtcheck 输出生成矩阵
    # gtcheck 输出格式: DC, sample1, sample2, discordance, concordance, nsites
    echo -e "sample1\\tsample2\\tconcordance\\tstatus" > identity_matrix.tsv
    touch identity_warnings.txt

    grep "^DC" gtcheck_output.txt | while read line; do
        S1=\$(echo "\$line" | awk '{print \$2}')
        S2=\$(echo "\$line" | awk '{print \$3}')
        DISCORD=\$(echo "\$line" | awk '{print \$4}')
        CONCORD=\$(echo "\$line" | awk '{print \$5}')
        NSITES=\$(echo "\$line" | awk '{print \$6}')
        
        if [ "\$NSITES" -gt 0 ]; then
            RATE=\$(echo "scale=4; \$CONCORD / \$NSITES" | bc)
        else
            RATE=0
        fi

        if (( \$(echo "\$RATE >= ${concordance_threshold}" | bc -l) )); then
            STATUS="SAME_SAMPLE"
            echo "WARNING: \$S1 and \$S2 appear to be the same sample (concordance: \$RATE)" >> identity_warnings.txt
        elif (( \$(echo "\$RATE >= 0.4" | bc -l) )); then
            STATUS="RELATED"
        else
            STATUS="DIFFERENT"
        fi

        echo -e "\$S1\\t\$S2\\t\$RATE\\t\$STATUS" >> identity_matrix.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
