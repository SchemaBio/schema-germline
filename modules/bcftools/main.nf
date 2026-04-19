// BCFTOOLS 模块
// 用途：变异检测、基因型计算、BAF 矩阵生成
// 包含的process：
//   BCFTOOLS_BAF_MATRIX - 基于 SNP 位点计算 BAF 分布矩阵

process BCFTOOLS_BAF_MATRIX {
    tag "BCFTOOLS_BAF_MATRIX on ${alignment.baseName}"
    label 'process_low'
    label 'automap'

    input:
        path alignment         // BAM 比对文件
        path alignment_index   // 比对文件索引 (.bai)
        path fasta             // 参考基因组 FASTA
        path fasta_fai         // 参考基因组索引 (.fai)
        path snp_positions     // SNP 位点文件 (VCF 或 BED 格式)
        val max_depth          // 最大测序深度 (默认 200)

    output:
        path "*.baf_matrix.tsv", emit: baf_matrix
        path "*.baf_matrix.json", emit: baf_json

    script:
    def depth = max_depth ?: 200
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam)$/, '')
    // 判断输入文件格式
    def regions_param = snp_positions.name.endsWith('.vcf') || snp_positions.name.endsWith('.vcf.gz')
        ? "-T ${snp_positions}"
        : "-R ${snp_positions}"
    """
    # 使用 bcftools mpileup 计算指定 SNP 位点的覆盖度
    # -a FORMAT/AD,FORMAT/DP: 输出等位基因深度和总深度
    # --max-depth: 限制最大深度以加快计算速度
    # -Q: base quality 阈值
    # -q: mapping quality 阈值
    bcftools mpileup \\
        -f ${fasta} \\
        ${regions_param} \\
        -a FORMAT/AD,FORMAT/DP \\
        --max-depth ${depth} \\
        -Q 20 \\
        -q 30 \\
        -Oz \\
        -o ${sample_id}.pileup.vcf.gz \\
        ${alignment}

    # 创建索引
    bcftools index -t ${sample_id}.pileup.vcf.gz

    # 从 VCF 提取 BAF 矩阵
    # BAF = B-allele depth / total depth
    # 对于二倍体，B-allele 通常指非参考等位基因
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD\\t%DP]\\n' ${sample_id}.pileup.vcf.gz | \\
    awk -F'\\t' '
    BEGIN {
        print "chrom\\tpos\\tref\\talt\\tad_ref\\tad_alt\\tdp\\tbaf"
    }
    {
        chrom = \$1
        pos = \$2
        ref = \$3
        alt = \$4
        # AD 格式: ref_depth,alt_depth (可能有多个 ALT)
        split(\$5, ad_arr, ",")
        ad_ref = ad_arr[1]
        ad_alt = ad_arr[2] + 0  # 如果没有 ALT 则为 0
        dp = \$6

        # 计算 BAF (B-Allele Frequency)
        if (dp > 0) {
            baf = ad_alt / dp
        } else {
            baf = 0
        }

        printf "%s\\t%s\\t%s\\t%s\\t%d\\t%d\\t%d\\t%.4f\\n", chrom, pos, ref, alt, ad_ref, ad_alt, dp, baf
    }
    ' > ${sample_id}.baf_matrix.tsv

    # 生成 JSON 格式输出
    echo "{" > ${sample_id}.baf_matrix.json
    echo "  \"sample_id\": \"${sample_id}\"," >> ${sample_id}.baf_matrix.json
    echo "  \"max_depth\": ${depth}," >> ${sample_id}.baf_matrix.json
    echo "  \"snp_count\": $(tail -n +2 ${sample_id}.baf_matrix.tsv | wc -l)," >> ${sample_id}.baf_matrix.json
    echo "  \"variants\": [" >> ${sample_id}.baf_matrix.json

    tail -n +2 ${sample_id}.baf_matrix.tsv | awk -F'\\t' '
    BEGIN { first = 1 }
    {
        if (!first) print ","
        first = 0
        printf "    {\"chrom\": \"%s\", \"pos\": %d, \"ref\": \"%s\", \"alt\": \"%s\", \"ad_ref\": %d, \"ad_alt\": %d, \"dp\": %d, \"baf\": %.4f}", \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8
    }
    ' >> ${sample_id}.baf_matrix.json

    echo "" >> ${sample_id}.baf_matrix.json
    echo "  ]" >> ${sample_id}.baf_matrix.json
    echo "}" >> ${sample_id}.baf_matrix.json

    # 清理临时文件
    rm -f ${sample_id}.pileup.vcf.gz ${sample_id}.pileup.vcf.gz.tbi
    """
}

