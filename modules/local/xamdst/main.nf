/*
 * Xamdst 覆盖度统计模块
 *
 * 功能：WES/WGS 样本的覆盖度统计
 * 工具：xamdst (bamdst 的现代替代品，支持 CRAM 和多线程)
 * 输出：覆盖度报告、区域统计、UTR/外显子/基因覆盖度、JSON 汇总
 */

process XAMDST {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(alignment), path(alignment_index)  // BAM/CRAM 及其索引
    path  target                      // 捕获区域 BED 文件 (外显子区域)
    path  target_fai                  // BED 索引 (可选)

    output:
    tuple val(meta), path("*.coverage.stat")     , emit: coverage_stat     // 覆盖度统计摘要 (文本)
    tuple val(meta), path("*.coverage.report.json"), emit: json             // JSON 格式汇总 (推荐用于下游)
    tuple val(meta), path("*.depth.stat")        , emit: depth_stat        // 深度分布统计
    tuple val(meta), path("*.chromosome.stat")   , emit: chromosome_stat   // 染色体覆盖度统计
    tuple val(meta), path("*.region.stat")       , emit: region_stat       // 区域覆盖度统计
    tuple val(meta), path("*.genelist.stat")     , emit: genelist_stat     // 基因列表覆盖度统计
    tuple val(meta), path("*.distribution.stat") , emit: distribution      // 覆盖度分布统计
    tuple val(meta), path("*.plot")              , emit: plot              // 覆盖度分布图 (可选)
    path "versions.yml"                          , emit: versions          // 软件版本

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def target_param = target ? "--target ${target}" : ''
    def cutoffdepth = params.xamdst_cutoffdepth ?: 10  // 默认 10X 覆盖度阈值
    """
    # 运行 xamdst 分析
    xamdst \
        --threads ${task.cpus} \
        --name ${prefix} \
        ${target_param} \
        --cutoffdepth ${cutoffdepth} \
        --bam ${alignment} \
        --outdir ./

    # 重命名输出文件
    if [ -f "coverage.report" ]; then mv "coverage.report" "${prefix}.coverage.stat"; fi
    if [ -f "depth distribution.data" ]; then mv "depth distribution.data" "${prefix}.depth.stat"; fi
    if [ -f "chromosome.stat" ]; then mv "chromosome.stat" "${prefix}.chromosome.stat"; fi
    if [ -f "region.stat" ]; then mv "region.stat" "${prefix}.region.stat"; fi
    if [ -f "genelist.stat" ]; then mv "genelist.stat" "${prefix}.genelist.stat"; fi
    if [ -f "depth distribution.info" ]; then mv "depth distribution.info" "${prefix}.distribution.stat"; fi
    if [ -f "coverage_distrib.png" ]; then mv "coverage_distrib.png" "${prefix}.plot"; fi

    # JSON 输出文件重命名
    if [ -f "coverage.report.json" ]; then mv "coverage.report.json" "${prefix}.coverage.report.json"; fi

    # 保留原始命名（兼容 xamdst 默认输出）
    if [ ! -f "${prefix}.coverage.stat" ] && [ -f "coverage.report" ]; then
        cp coverage.report ${prefix}.coverage.stat
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xamdst: \$(xamdst --version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}
