/*
 * SVDB - Structural Variation Database Builder & Query
 *
 * 功能：结构变异数据库构建、合并、查询和注释
 * 工具：SVDB
 *
 * 说明：
 *   - 支持 VCF 文件合并
 *   - 支持 SV 数据库构建
 *   - 支持 CNV 注释（频率、来源等）
 */

process SVDB_QUERY {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)  // 输入 VCF (CNV/SV)
    path  db_vcf                                // 数据库 VCF (可选)
    path  db_vcf_idx                            // 数据库索引 (可选)
    val   db_type                               // 数据库类型: 'gnomad', 'sweegen', 'custom'
    val   in_occ                                // 数据库中的计数标签
    val   in_frq                                // 数据库中的频率标签
    val   out_occ                               // 输出的计数标签名
    val   out_frq                               // 输出的频率标签名
    val   bnd_distance                          // 断点距离阈值
    val   overlap                               // 重叠比例阈值

    output:
    tuple val(meta), path("*.svdb.vcf.gz"), path("*.svdb.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_cmd = db_vcf.name != 'NO_FILE' ? "--db ${db_vcf}" : ''
    def in_occ_cmd = in_occ ? "--in_occ ${in_occ}" : ''
    def in_frq_cmd = in_frq ? "--in_frq ${in_frq}" : ''
    def out_occ_cmd = out_occ ? "--out_occ ${out_occ}" : ''
    def out_frq_cmd = out_frq ? "--out_frq ${out_frq}" : ''
    def bnd_cmd = bnd_distance ? "--bnd_distance ${bnd_distance}" : '--bnd_distance 2500'
    def overlap_cmd = overlap ? "--overlap ${overlap}" : '--overlap 0.8'
    """
    svdb --query \\
        --query_vcf ${vcf} \\
        ${db_cmd} \\
        ${in_occ_cmd} \\
        ${in_frq_cmd} \\
        ${out_occ_cmd} \\
        ${out_frq_cmd} \\
        ${bnd_cmd} \\
        ${overlap_cmd} \\
        --prefix ${prefix}.svdb \\
        ${args}

    # 压缩并索引
    bgzip -f ${prefix}.svdb.vcf
    tabix -f -p vcf ${prefix}.svdb.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$(svdb --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * SVDB_MERGE - 合并多个 VCF 文件
 *
 * 功能：合并来自不同 SV/CNV calling 工具的 VCF
 */
process SVDB_MERGE {
    label 'process_medium'

    input:
    path vcfs                     // 多个 VCF 文件
    val  bnd_distance             // 断点距离阈值
    val  overlap                  // 重叠比例阈值
    val  priority                 // 优先级列表 (逗号分隔)
    val  pass_only                // 只合并 PASS 变体

    output:
    path "merged.vcf.gz", emit: vcf
    path "merged.vcf.gz.tbi", emit: tbi
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def bnd_cmd = bnd_distance ? "--bnd_distance ${bnd_distance}" : '--bnd_distance 10000'
    def overlap_cmd = overlap ? "--overlap ${overlap}" : '--overlap 0.6'
    def prio_cmd = priority ? "--priority ${priority}" : ''
    def pass_cmd = pass_only ? '--pass_only' : ''

    // 将 VCF 列表转换为空格分隔的参数
    def vcf_list = vcfs.collect { "--vcf ${it}" }.join(' ')

    """
    svdb --merge \\
        ${vcf_list} \\
        ${bnd_cmd} \\
        ${overlap_cmd} \\
        ${prio_cmd} \\
        ${pass_cmd} \\
        ${args} | bgzip -c > merged.vcf.gz

    tabix -f -p vcf merged.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$(svdb --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * SVDB_BUILD - 构建 SV 数据库
 *
 * 功能：从多个 VCF 文件构建 SQLite 数据库
 */
process SVDB_BUILD {
    label 'process_high'

    input:
    path vcfs                      // 多个 VCF 文件
    path folder                    // 或文件夹 (二选一)
    val   prefix                   // 输出前缀
    val   bnd_distance             // 断点距离阈值
    val   overlap                  // 重叠比例阈值
    val   use_dbscan               // 使用 DBSCAN 聚类
    val   epsilon                  // DBSCAN epsilon
    val   min_pts                  // DBSCAN min_pts

    output:
    path "*.db", emit: database
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix_cmd = prefix ? "--prefix ${prefix}" : ''
    def bnd_cmd = bnd_distance ? "--bnd_distance ${bnd_distance}" : ''
    def overlap_cmd = overlap ? "--overlap ${overlap}" : ''
    def dbscan_cmd = use_dbscan ? '--DBSCAN' : ''
    def epsilon_cmd = epsilon ? "--epsilon ${epsilon}" : '--epsilon 500'
    def min_pts_cmd = min_pts ? "--min_pts ${min_pts}" : '--min_pts 2'

    // VCF 输入
    if (vcfs) {
        def vcf_list = vcfs.collect { "${it}" }.join(' ')
        """
        svdb --build \\
            ${vcf_list} \\
            ${prefix_cmd} \\
            ${bnd_cmd} \\
            ${overlap_cmd} \\
            ${dbscan_cmd} \\
            ${epsilon_cmd} \\
            ${min_pts_cmd} \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$(svdb --version 2>&1 || echo "unknown")
        END_VERSIONS
        """
    } else {
        """
        svdb --build \\
            --folder ${folder} \\
            ${prefix_cmd} \\
            ${bnd_cmd} \\
            ${overlap_cmd} \\
            ${dbscan_cmd} \\
            ${epsilon_cmd} \\
            ${min_pts_cmd} \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            svdb: \$(svdb --version 2>&1 || echo "unknown")
        END_VERSIONS
        """
    }
}

/*
 * SVDB_EXPORT - 导出数据库变体
 *
 * 功能：从数据库导出变体为 VCF
 */
process SVDB_EXPORT {
    label 'process_medium'

    input:
    path  database                 // SVDB 数据库文件 (.db)
    val   bnd_distance             // 断点距离阈值
    val   overlap                  // 重叠比例阈值
    val   use_dbscan               // 使用 DBSCAN
    val   prefix                   // 输出前缀

    output:
    path "*.vcf.gz", emit: vcf
    path "*.vcf.gz.tbi", emit: tbi
    path "versions.yml"    , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix_cmd = prefix ? "--prefix ${prefix}" : ''
    def bnd_cmd = bnd_distance ? "--bnd_distance ${bnd_distance}" : '--bnd_distance 2500'
    def overlap_cmd = overlap ? "--overlap ${overlap}" : '--overlap 0.8'
    def dbscan_cmd = use_dbscan ? '--DBSCAN' : ''

    """
    svdb --export \\
        --db ${database} \\
        ${prefix_cmd} \\
        ${bnd_cmd} \\
        ${overlap_cmd} \\
        ${dbscan_cmd} \\
        ${args} | bgzip -c > ${prefix ?: 'exported'}.vcf.gz

    tabix -f -p vcf ${prefix ?: 'exported'}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$(svdb --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}
