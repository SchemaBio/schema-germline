// GATK 变异检测模块
// 用途：标记重复、合并比对、质控统计、变异检测
// 包含的process：
//   MARKDUPLICATES - 标记PCR重复
//   COLLECTQCMETRICS - 收集质控指标
//   MUTECT2_MT - 线粒体变异检测 (线粒体模式)

process MARKDUPLICATES {
    tag "MARKDUPLICATES on $sample_id"
    label 'gatk'
    label 'process_medium'

    input:
        path alignment      // BAM 文件 (文件名包含 sample_id)
        path alignment_index // BAM 索引文件
        val fasta           // 参考基因组路径

    output:
        path "*.marked.bam", emit: alignment
        path("*.marked.metrics.txt"), emit: metrics

    script:
    // 从文件名提取 sample_id
    def sample_id = alignment.baseName.replaceAll(/\.bam$/, '')
    """
    gatk MarkDuplicates \\
        -I ${alignment} \\
        -O ${sample_id}.marked.bam \\
        -M ${sample_id}.marked.metrics.txt \\
        --CREATE_INDEX false \\
        --REFERENCE_SEQUENCE ${fasta}
    """
}

process COLLECTQCMETRICS {
    tag "COLLECTQCMETRICS on ${alignment.baseName}"
    label 'gatk'
    label 'process_medium'

    input:
        path alignment
        path alignment_index
        path fasta
        path target_bed

    output:
        path "*.metrics", emit: metrics
        path "*.pdf", emit: pdf, optional: true

    script:
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam)$/, '')
    """
    gatk CollectMultipleMetrics \\
        -I ${alignment} \\
        -O ${sample_id} \\
        -R ${fasta} \\
        --INTERVALS ${target_bed} \\
        --PROGRAM CollectAlignmentSummaryMetrics \\
        --PROGRAM CollectInsertSizeMetrics \\
        --PROGRAM QualityScoreDistribution \\
        --PROGRAM CollectGcBiasMetrics \\
        --PROGRAM MeanQualityByCycle \\
        --PROGRAM CollectBaseDistributionByCycle
    """
}

process MUTECT2_MT {
    tag "MUTECT2_MT on ${alignment.baseName}"
    label 'gatk'
    label 'process_medium'

    input:
        path alignment           // BAM 比对文件
        path alignment_index     // 比对文件索引 (.bai)
        path fasta               // 参考基因组 FASTA (包含线粒体序列)
        path fasta_fai           // 参考基因组索引
        path fasta_dict          // 参考基因组字典
        val genome_assembly      // 基因组版本: 'GRCh37' 或 'GRCh38'

    output:
        path "*.mt.vcf.gz", emit: vcf
        path "*.mt.vcf.gz.tbi", emit: vcf_tbi
        path "*.mt.vcf.gz.stats", emit: stats

    script:
    def sample_id = alignment.baseName.replaceAll(/\.(marked|bam)$/, '')
    // 线粒体染色体名称映射
    // GRCh37/GRCh38 可能使用不同的命名: MT, chrM, chrMT, M
    """
    # 检测线粒体染色体名称
    # 优先顺序: chrM > MT > chrMT > M
    MT_CHR=""
    for chr in "chrM" "MT" "chrMT" "M"; do
        if samtools view -H ${alignment} 2>/dev/null | grep -q "@SQ.*SN:\${chr}"; then
            MT_CHR="\${chr}"
            break
        fi
    done

    if [ -z "\${MT_CHR}" ]; then
        echo "Warning: No mitochondrial chromosome found in reference" >&2
        # 创建空 VCF
        echo "##fileformat=VCFv4.2" > ${sample_id}.mt.vcf
        echo "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" >> ${sample_id}.mt.vcf
        echo "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" >> ${sample_id}.mt.vcf
        echo "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${sample_id}" >> ${sample_id}.mt.vcf
        bgzip -c ${sample_id}.mt.vcf > ${sample_id}.mt.vcf.gz
        tabix -p vcf ${sample_id}.mt.vcf.gz
        exit 0
    fi

    echo "Detected mitochondrial chromosome: \${MT_CHR}"

    # 运行 Mutect2 线粒体模式
    # --mitochondria-mode: 自动调整参数适配 mtDNA
    # --max-reads-per-alignment-start 0: 不限制 reads 数量
    # --max-mnp-distance 0: 不合并相邻变异
    gatk Mutect2 \\
        -R ${fasta} \\
        -I ${alignment} \\
        -O ${sample_id}.mt.unfiltered.vcf.gz \\
        --mitochondria-mode \\
        --max-reads-per-alignment-start 0 \\
        --max-mnp-distance 0 \\
        -L \${MT_CHR}

    # 过滤变异 (线粒体模式使用简单过滤)
    gatk FilterMutectCalls \\
        -V ${sample_id}.mt.unfiltered.vcf.gz \\
        -R ${fasta} \\
        -O ${sample_id}.mt.vcf.gz

    # 清理中间文件
    rm -f ${sample_id}.mt.unfiltered.vcf.gz ${sample_id}.mt.unfiltered.vcf.gz.tbi ${sample_id}.mt.unfiltered.vcf.gz.stats
    """
}