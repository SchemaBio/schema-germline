// Sambamba 模块
// 用途：BAM 文件的重复标记
// 包含的process：
//   SAMBAMBA_MARKDUP - 标记或移除 PCR 重复 reads
//
// 与 GATK MarkDuplicates 相比的优势：
//   - 更快的处理速度 (多线程优化)
//   - 更低的内存占用
//   - 自动创建索引

process SAMBAMBA_MARKDUP {
    tag "SAMBAMBA_MARKDUP on ${alignment.baseName}"
    label 'mapping'
    label 'process_medium'

    input:
        path alignment          // BAM 文件
        path alignment_index    // BAM 索引文件 (.bai)

    output:
        path "*.marked.bam", emit: alignment
        path "*.marked.bam.bai", emit: index
        path "*.marked.metrics.txt", emit: metrics

    script:
    def threads = Math.min(task.cpus as int, 8)
    def sample_id = alignment.baseName.replaceAll(/\.bam$/, '')
    def remove_dup = params.sambamba_remove_duplicates ?: false
    def remove_flag = remove_dup ? '--remove-duplicates' : ''
    """
    sambamba markdup \\
        -t ${threads} \\
        ${remove_flag} \\
        --tmpdir=. \\
        --overflow-list-size=500000 \\
        --hash-table-size=500000 \\
        ${alignment} ${sample_id}.marked.bam

    # 创建索引
    sambamba index -t ${threads} ${sample_id}.marked.bam

    # 生成 metrics 文件 (类似 GATK 格式)
    # sambamba markdup 不直接输出 metrics，使用 samtools stats 计算
    samtools stats ${sample_id}.marked.bam > ${sample_id}.marked.stats.txt

    # 提取重复相关统计并转换为类似 GATK metrics 格式
    TOTAL_READS=\$(grep "^SN" ${sample_id}.marked.stats.txt | grep "raw total sequences" | cut -d: -f2 | tr -d ' ')
    DUP_READS=\$(grep "^SN" ${sample_id}.marked.stats.txt | grep "reads marked as duplicates" | cut -d: -f2 | tr -d ' ')
    DUP_PERCENT=\$(echo "scale=2; \${DUP_READS} * 100 / \${TOTAL_READS}" | bc 2>/dev/null || echo "0")

    cat <<EOF > ${sample_id}.marked.metrics.txt
## METRICS CLASS        Picard.metrics.MarkDuplicates
READ_PAIRS_EXAMINED	${TOTAL_READS}
UNPAIRED_READS_EXAMINED	0
UNPAIRED_READ_DUPLICATES	0
READ_PAIR_DUPLICATES	${DUP_READS}
READ_PAIR_DUPLICATE_PERCENT	${DUP_PERCENT}
PERCENT_DUPLICATION	${DUP_PERCENT}
EOF
    """
}