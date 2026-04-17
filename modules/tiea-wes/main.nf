// TIEA-WES 模块
// 用途：移动元件插入 (MEI) 检测
// 工具：TIEA-WES (Transposon Insertion Event Analyzer for WES)
// 输出：VCF 格式的 TE 插入结果 (Alu/LINE-1/SVA/HERV)
//
// 说明：
//   - 检测 Alu、LINE-1、SVA、HERV 等转座子插入事件
//   - 基于 softclip 信号的 bidirectional 检测 (5' 和 3')
//   - 支持 BAM 输入格式
//   - Docker 镜像已包含 TE 基因组参考 (hg38reps.fa) 和配置

// ============================================================================
// TIEA_WES - MEI 检测
// ============================================================================
process TIEA_WES {
    tag "TIEA_WES on ${alignment.baseName}"
    label 'process_low'
    label 'tiea_wes'

    input:
        path alignment           // BAM 比对文件
        path alignment_index     // 比对文件索引
        path fasta               // 参考基因组 FASTA
        val sample_id            // 样本标识符
        val min_support          // 最小断点支持 reads 数 (默认 10)
        val min_softclip         // 最小 softclip 长度 (默认 36)
        val min_mapq             // 最小 MAPQ (默认 20)
        val cluster_window       // 断点聚类窗口 bp (默认 10)
        val threads              // BWA 线程数 (默认 4)

    output:
        path "${sample_id}.te.result.vcf.gz", emit: vcf
        path "${sample_id}.te.result.vcf.gz.tbi", emit: vcf_tbi

    script:
    def support = min_support ?: 10
    def softclip = min_softclip ?: 36
    def mapq = min_mapq ?: 20
    def window = cluster_window ?: 10
    def num_threads = threads ?: 4
    def ref_cmd = fasta.name != 'NO_FILE' ? "-r ${fasta}" : ''
    """
    # 运行 TIEA-WES
    # Docker 镜像中已配置好 TE reference 和 config.ini
    python /app/TIEA-WES.py \\
        -p ${sample_id} \\
        -i ${alignment} \\
        -o . \\
        ${ref_cmd} \\
        -s ${support} \\
        -l ${softclip} \\
        -q ${mapq} \\
        -w ${window} \\
        --threads ${num_threads}

    # 压缩 VCF 并创建索引
    bgzip -c "${sample_id}.te.result.vcf" > "${sample_id}.te.result.vcf.gz"
    tabix -p vcf "${sample_id}.te.result.vcf.gz"

    # 清理原始 VCF
    rm -f "${sample_id}.te.result.vcf"
    """
}