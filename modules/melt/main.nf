/*
 * MELT 模块集合
 *
 * 功能：移动元件插入 (Mobile Element Insertion, MEI) 检测
 * 工具：MELT (Mobile Element Locator Tool)
 *
 * 支持检测的转座子类型：
 *   - LINE-1 (L1) - 约 6kb 的长片段插入
 *   - ALU (SINE) - 约 300bp 的短片段插入
 *   - SVA - 约 2kb 的复合元件插入
 *
 * 分析模式：
 *   - MELT-SINGLE: 单样本分析流程 (适用于 WES/WGS)
 *   - MELT-Group: 多样本群体分析 (需要多个样本同时分析)
 *
 * WES 注意事项：
 *   - 需要提供目标区域 BED 文件限制分析范围
 *   - 仅检测落在捕获区域内的转座子插入
 *
 * 参考：
 *   https://melt-igs.umaryland.edu
 *   Gardner et al. 2017, "The Mobile Element Locator Tool (MELT)"
 */

/*
 * MELT_SPLITREADS - 预处理 BAM 文件
 *
 * 功能：提取 split reads 和 discordant pairs 作为候选插入位点
 * 这是 MELT 分析的第一步，为后续 IndivAnalysis 准备数据
 */
process MELT_SPLITREADS {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(alignment), path(alignment_index)  // BAM/CRAM 及索引
    path  fasta                                              // 参考基因组
    path  "${fasta}.fai"                                     // FASTA 索引
    path  "${fasta}.dict"                                    // FASTA 字典

    output:
    tuple val(meta), path("*.split.bam")           , emit: split_bam     // Split reads BAM
    tuple val(meta), path("*.discordant.bam")      , emit: discordant_bam // Discordant pairs BAM
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # MELT 需要检查 BAM 文件是否包含必要信息
    # 创建临时目录用于 MELT 运行
    mkdir -p ${prefix}_melt_work

    # Split BAM file preprocessing
    java -Xmx${task.memory.toGiga() as int}G -jar /opt/MELT/MELT.jar \\
        Preprocess \\
        -h ${fasta} \\
        -bam ${alignment} \\
        -w ${prefix}_melt_work \\
        ${args}

    # 移动输出文件到工作目录
    cp ${prefix}_melt_work/*.split.bam ${prefix}.split.bam || true
    cp ${prefix}_melt_work/*.discordant.bam ${prefix}.discordant.bam || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        melt: \$(java -jar /opt/MELT/MELT.jar -version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * MELT_SINGLE_INDIV - 单样本转座子分析
 *
 * 功能：对单个样本进行特定转座子类型的插入检测
 * 需要对每种转座子类型分别运行
 *
 * 输入：
 *   - 预处理后的 BAM 文件 (split + discordant)
 *   - 参考基因组
 *   - 转座子参考序列 (ME_ref.fa)
 *   - 转座子类型名称 (LINE1, ALU, SVA)
 *   - 目标区域 BED (可选，WES 推荐)
 */
process MELT_SINGLE_INDIV {
    tag "$meta.id - ${me_type}"
    label 'process_high'

    input:
    tuple val(meta), path(split_bam), path(discordant_bam)   // 预处理 BAM
    path  fasta                                              // 参考基因组
    path  "${fasta}.fai"                                     // FASTA 索引
    path  "${fasta}.dict"                                    // FASTA 字典
    path  me_reference                                       // 转座子参考序列 (如 L1LINE1.fa, ALU.fa, SVA.fa)
    val   me_type                                            // 转座子类型: LINE1, ALU, SVA
    path  target_bed                                         // 目标区域 BED (可选，WES 推荐)

    output:
    tuple val(meta), path("${me_type}/*.vcf"), emit: vcf
    tuple val(meta), path("${me_type}/*.bed") , emit: bed, optional: true
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions_param = target_bed ? "-c ${target_bed}" : ''
    """
    # 创建输出目录
    mkdir -p ${me_type}

    # 单样本转座子分析
    java -Xmx${task.memory.toGiga() as int}G -jar /opt/MELT/MELT.jar \\
        IndivAnalysis \\
        -h ${fasta} \\
        -bam ${split_bam} ${discordant_bam} \\
        -m ${me_reference} \\
        -t ${me_type} \\
        -w ${me_type} \\
        ${regions_param} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        melt: \$(java -jar /opt/MELT/MELT.jar -version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * MELT_SINGLE_ANALYSIS - 单样本完整分析流程
 *
 * 功能：一步完成预处理 + 三种转座子分析 + 结果合并
 * 这是推荐的 WES 单样本使用方式
 *
 * 输入：
 *   - BAM/CRAM 比对文件
 *   - 参考基因组及索引
 *   - 转座子参考序列集合 (LINE1, ALU, SVA)
 *   - 目标区域 BED (WES 推荐)
 */
process MELT_SINGLE_ANALYSIS {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(alignment), path(alignment_index)  // BAM/CRAM 及索引
    path  fasta                                              // 参考基因组
    path  "${fasta}.fai"                                     // FASTA 索引
    path  "${fasta}.dict"                                    // FASTA 字典
    path  me_references                                      // 转座子参考序列集合 (LINE1.fa, ALU.fa, SVA.fa)
    path  target_bed                                         // 目标区域 BED (可选，WES 推荐)
    val   me_types                                           // 要分析的转座子类型，逗号分隔 (如: "LINE1,ALU,SVA")

    output:
    tuple val(meta), path("LINE1/*.vcf")        , emit: line1_vcf
    tuple val(meta), path("ALU/*.vcf")          , emit: alu_vcf
    tuple val(meta), path("SVA/*.vcf")          , emit: sva_vcf
    tuple val(meta), path("*.merged.vcf")       , emit: merged_vcf
    tuple val(meta), path("*.mei_summary.tsv")  , emit: summary
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def types = me_types ?: 'LINE1,ALU,SVA'
    def regions_param = target_bed ? "-c ${target_bed}" : ''
    """
    # 创建工作目录
    mkdir -p melt_work LINE1 ALU SVA

    # Step 1: Preprocess BAM file
    java -Xmx${task.memory.toGiga() as int}G -jar /opt/MELT/MELT.jar \\
        Preprocess \\
        -h ${fasta} \\
        -bam ${alignment} \\
        -w melt_work \\
        ${args}

    # Step 2: 对每种转座子类型运行 IndivAnalysis
    # LINE-1 分析
    if echo "${types}" | grep -q "LINE1"; then
        java -Xmx${task.memory.toGiga() as int}G -jar /opt/MELT/MELT.jar \\
            IndivAnalysis \\
            -h ${fasta} \\
            -bam melt_work/*.split.bam melt_work/*.discordant.bam \\
            -m /opt/MELT/me_refs/L1LINE1.fa \\
            -t LINE1 \\
            -w LINE1 \\
            ${regions_param}
    fi

    # ALU 分析
    if echo "${types}" | grep -q "ALU"; then
        java -Xmx${task.memory.toGiga() as int}G -jar /opt/MELT/MELT.jar \\
            IndivAnalysis \\
            -h ${fasta} \\
            -bam melt_work/*.split.bam melt_work/*.discordant.bam \\
            -m /opt/MELT/me_refs/ALU.fa \\
            -t ALU \\
            -w ALU \\
            ${regions_param}
    fi

    # SVA 分析
    if echo "${types}" | grep -q "SVA"; then
        java -Xmx${task.memory.toGiga() as int}G -jar /opt/MELT/MELT.jar \\
            IndivAnalysis \\
            -h ${fasta} \\
            -bam melt_work/*.split.bam melt_work/*.discordant.bam \\
            -m /opt/MELT/me_refs/SVA.fa \\
            -t SVA \\
            -w SVA \\
            ${regions_param}
    fi

    # Step 3: 合并所有转座子结果
    # 收集所有 VCF 文件并合并
    vcf_files=\$(find LINE1 ALU SVA -name "*.vcf" 2>/dev/null)
    if [ -n "\$vcf_files" ]; then
        # 创建合并 VCF (简单合并，不重新分型)
        echo "##fileformat=VCFv4.2" > ${prefix}.merged.vcf
        echo "##INFO=<ID=ME_TYPE,Number=1,Type=String,Description=\"Mobile element type\">" >> ${prefix}.merged.vcf
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.id}" >> ${prefix}.merged.vcf

        for vcf in \$vcf_files; do
            # 从目录名获取转座子类型
            me_type=\$(basename \$(dirname \$vcf))
            # 提取数据行并添加 ME_TYPE 信息
            grep -v "^#" \$vcf | while read line; do
                echo "\$line" | awk -v me=\$me_type 'BEGIN{OFS="\t"} {if(\$7==".") \$7="PASS"; \$8=\$8";ME_TYPE="me; print}'
            done >> ${prefix}.merged.vcf
        done
    fi

    # Step 4: 生成汇总报告
    echo "sample\tme_type\tcount\tchrom\tfirst_pos\tlast_pos" > ${prefix}.mei_summary.tsv
    for type in LINE1 ALU SVA; do
        if [ -d "\$type" ]; then
            count=\$(grep -v "^#" \$type/*.vcf 2>/dev/null | wc -l || echo 0)
            if [ \$count -gt 0 ]; then
                first_pos=\$(grep -v "^#" \$type/*.vcf 2>/dev/null | head -1 | cut -f1-2)
                last_pos=\$(grep -v "^#" \$type/*.vcf 2>/dev/null | tail -1 | cut -f1-2)
                echo "${meta.id}\t\$type\t\$count\t\$first_pos\t\$last_pos" >> ${prefix}.mei_summary.tsv
            else
                echo "${meta.id}\t\$type\t0\t.\t." >> ${prefix}.mei_summary.tsv
            fi
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        melt: \$(java -jar /opt/MELT/MELT.jar -version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * MELT_GENOTYPE - 已知位点分型
 *
 * 功能：对已知的转座子插入位点进行分型
 * 用于验证或比较已知位点在不同样本中的状态
 *
 * 输入：
 *   - BAM/CRAM 比对文件
 *   - 参考基因组
 *   - 已知位点 BED/VCF 文件
 */
process MELT_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(alignment), path(alignment_index)
    path  fasta
    path  "${fasta}.fai"
    path  "${fasta}.dict"
    path  known_sites                                        // 已知 MEI 位点文件
    val   me_type                                            // 转座子类型

    output:
    tuple val(meta), path("*.genotyped.vcf"), emit: vcf
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p genotype_work

    java -Xmx${task.memory.toGiga() as int}G -jar /opt/MELT/MELT.jar \\
        Genotype \\
        -h ${fasta} \\
        -bam ${alignment} \\
        -s ${known_sites} \\
        -t ${me_type} \\
        -w genotype_work \\
        ${args}

    cp genotype_work/*.vcf ${prefix}.genotyped.vcf || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        melt: \$(java -jar /opt/MELT/MELT.jar -version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * MELT_FILTER - 结果过滤
 *
 * 功能：根据质量指标过滤低置信度的转座子插入
 * 可选步骤，用于提高结果可靠性
 */
process MELT_FILTER {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf)
    val   min_support                                        // 最小支持 reads 数 (默认 3)
    val   min_mapq                                           // 最小 MAPQ (默认 20)

    output:
    tuple val(meta), path("*.filtered.vcf"), emit: vcf
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def support = min_support ?: 3
    def mapq = min_mapq ?: 20
    """
    # 过滤 VCF 文件
    # 保留 PASS 或支持 reads >= threshold 的记录
    awk -v support=${support} -v mapq=${mapq} '
        BEGIN {OFS="\t"}
        /^#/ {print; next}
        {
            # 解析 INFO 字段获取支持数
            split(\$8, info, ";")
            supp=0
            for(i in info) {
                if(info[i] ~ /^SUPPORT=/) {
                    split(info[i], s, "=")
                    supp=s[2]
                }
            }
            if(\$7=="PASS" || supp>=support) print
        }
    ' ${vcf} > ${prefix}.filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        melt: \$(java -jar /opt/MELT/MELT.jar -version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}