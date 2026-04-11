// VEP (Ensembl Variant Effect Predictor) 模块
// 用途：对 VCF 文件进行变异功能注释
// 工具：Ensembl VEP release 115.2
// 输出：注释后的 VCF 文件，包含基因、转录本、蛋白等注释信息
//
// 主要参数说明：
//   - cache: 使用本地缓存加速注释
//   - offline: 离线模式，不查询在线数据库
//   - pick: 选择每个变异最显著的注释结果
//   - everything: 输出所有可用注释字段
//   - fork: 多线程并行处理

process VEP_ANNOTATE {
    tag "VEP on ${vcf.baseName}"
    label 'process_medium'
    label 'vep'
    publishDir "${params.output}/05.Annotations/SNP_InDel", mode: 'copy'

    input:
        path vcf               // 输入 VCF 文件
        path vcf_tbi           // VCF 索引文件 (.tbi)
        path fasta             // 参考基因组 FASTA
        path fasta_fai         // 参考基因组索引 (.fai)
        val genome_assembly    // 基因组版本: GRCh37 或 GRCh38
        val use_pick           // 是否选择最显著注释 (默认 true)
        val use_refseq_only    // 是否仅使用 RefSeq 转录本 (默认 false)
        val cache_dir          // VEP 缓存目录 (可选，默认使用内置缓存)
        val extra_args         // 额外参数 (可选)

    output:
        path "*.vep.vcf.gz", emit: vep_vcf
        path "*.vep.vcf.gz.tbi", emit: vep_vcf_tbi
        path "*.vep.summary.txt", emit: summary

    script:
    def assembly = genome_assembly ?: 'GRCh38'
    def pick = use_pick ? '--pick' : ''
    def refseq = use_refseq_only ? '--refseq' : ''
    def cache = cache_dir ? "--cache --dir_cache ${cache_dir}" : '--cache'
    def extra = extra_args ?: ''
    def sample_id = vcf.baseName.replaceAll(/\.(vcf|vcf\.gz)$/, '')
    def threads = Math.min(task.cpus as int, 8)  // VEP fork 最多 8 线程较稳定
    """
    # 设置 VEP 缓存目录
    # 如果未指定外部缓存目录，使用容器内置缓存
    if [ "${cache_dir}" != "null" ] && [ -d "${cache_dir}" ]; then
        VEP_CACHE="${cache_dir}"
    else
        # 使用容器内置缓存路径
        VEP_CACHE="/opt/vep/.vep"
    fi

    # 运行 VEP 注释
    vep \\
        --input_file ${vcf} \\
        --output_file ${sample_id}.vep.vcf \\
        --format vcf \\
        --vcf \\
        --assembly ${assembly} \\
        ${cache} \\
        --offline \\
        --dir_cache \${VEP_CACHE} \\
        --fasta ${fasta} \\
        ${pick} \\
        ${refseq} \\
        --everything \\
        --force_overwrite \\
        --no_stats \\
        --fork ${threads} \\
        ${extra}

    # 压缩 VCF 文件并创建索引
    bgzip -c ${sample_id}.vep.vcf > ${sample_id}.vep.vcf.gz
    tabix -p vcf ${sample_id}.vep.vcf.gz

    # 生成注释摘要
    echo "VEP Annotation Summary for ${sample_id}" > ${sample_id}.vep.summary.txt
    echo "Input VCF: ${vcf}" >> ${sample_id}.vep.summary.txt
    echo "Output VCF: ${sample_id}.vep.vcf.gz" >> ${sample_id}.vep.summary.txt
    echo "Assembly: ${assembly}" >> ${sample_id}.vep.summary.txt
    echo "Total variants annotated:" >> ${sample_id}.vep.summary.txt
    zcat ${sample_id}.vep.vcf.gz | grep -v "^#" | wc -l >> ${sample_id}.vep.summary.txt
    """
}


// VEP_MT - 线粒体变异注释
// 用途：对线粒体 VCF 文件进行特化注释
// 特点：使用 VEP 的 mitochondria 模式，支持线粒体特异性注释字段
// 输出：注释后的线粒体 VCF，包含 mtDNA 特异的功能预测
//
// 线粒体特化参数：
//   - mitochondria: 启用线粒体模式
//   - variant_class: 输出变异类型分类
//   - protein: 输出蛋白级别变化
//   - polyphen/sift: 功能预测评分

process VEP_MT {
    tag "VEP_MT on ${vcf.baseName}"
    label 'process_low'
    label 'vep'
    publishDir "${params.output}/05.Annotations/Mitochondria", mode: 'copy'

    input:
        path vcf               // 输入线粒体 VCF 文件
        path vcf_tbi           // VCF 索引文件 (.tbi)
        path fasta             // 参考基因组 FASTA
        path fasta_fai         // 参考基因组索引 (.fai)
        val genome_assembly    // 基因组版本: GRCh37 或 GRCh38
        val cache_dir          // VEP 缓存目录 (可选)
        val extra_args         // 额外参数 (可选)

    output:
        path "*.vep.mt.vcf.gz", emit: vep_vcf
        path "*.vep.mt.vcf.gz.tbi", emit: vep_vcf_tbi
        path "*.vep.mt.summary.txt", emit: summary

    script:
    def assembly = genome_assembly ?: 'GRCh38'
    def cache = cache_dir ? "--cache --dir_cache ${cache_dir}" : '--cache'
    def extra = extra_args ?: ''
    def sample_id = vcf.baseName.replaceAll(/\.(vcf|vcf\.gz)$/, '')
    """
    # 设置 VEP 缓存目录
    if [ "${cache_dir}" != "null" ] && [ -d "${cache_dir}" ]; then
        VEP_CACHE="${cache_dir}"
    else
        VEP_CACHE="/opt/vep/.vep"
    fi

    # 运行 VEP 线粒体注释
    # --mitochondria: 启用线粒体特化模式
    # 该模式会输出线粒体特异性注释字段，如：
    #   - MITOTIP: tRNA 功能预测
    #   - HmtVar: 线粒体变异数据库注释
    vep \\
        --input_file ${vcf} \\
        --output_file ${sample_id}.vep.mt.vcf \\
        --format vcf \\
        --vcf \\
        --assembly ${assembly} \\
        ${cache} \\
        --offline \\
        --dir_cache \${VEP_CACHE} \\
        --fasta ${fasta} \\
        --mitochondria \\
        --variant_class \\
        --protein \\
        --polyphen \\
        --sift \\
        --force_overwrite \\
        --no_stats \\
        ${extra}

    # 压缩 VCF 文件并创建索引
    bgzip -c ${sample_id}.vep.mt.vcf > ${sample_id}.vep.mt.vcf.gz
    tabix -p vcf ${sample_id}.vep.mt.vcf.gz

    # 生成线粒体注释摘要
    echo "VEP Mitochondria Annotation Summary for ${sample_id}" > ${sample_id}.vep.mt.summary.txt
    echo "Input VCF: ${vcf}" >> ${sample_id}.vep.mt.summary.txt
    echo "Output VCF: ${sample_id}.vep.mt.vcf.gz" >> ${sample_id}.vep.mt.summary.txt
    echo "Assembly: ${assembly}" >> ${sample_id}.vep.mt.summary.txt
    echo "Mitochondria variants annotated:" >> ${sample_id}.vep.mt.summary.txt
    zcat ${sample_id}.vep.mt.vcf.gz | grep -v "^#" | wc -l >> ${sample_id}.vep.mt.summary.txt
    echo "" >> ${sample_id}.vep.mt.summary.txt
    echo "Heteroplasmy variants (ALT frequency):" >> ${sample_id}.vep.mt.summary.txt
    zcat ${sample_id}.vep.mt.vcf.gz | grep -v "^#" | head -20 >> ${sample_id}.vep.mt.summary.txt
    """
}


// VEP_MEI - 移动元件插入 (MEI) 注释
// 用途：对 TIEA-WES 产生的 MEI VCF 进行功能注释
// 输入：TIEA-WES 输出的 TE 插入 VCF (Alu/LINE-1/SVA/HERV)
// 输出：注释后的 MEI VCF，包含插入位点附近的基因影响信息
//
// MEI 注释特点：
//   - MEI 是大片段插入事件，VEP 会标注 insertion 类型
//   - 注释关注插入位点对附近基因的影响（如 gene_upstream, gene_downstream）
//   - INFO 字段中的 RN 标识转座子类型 (ALU/L1/SVA/HERV)
//   - 使用 --distance 调整上下游基因距离阈值

process VEP_MEI {
    tag "VEP_MEI on ${vcf.baseName}"
    label 'process_low'
    label 'vep'
    publishDir "${params.output}/05.Annotations/MEI", mode: 'copy'

    input:
        path vcf               // 输入 MEI VCF 文件 (TIEA-WES 输出)
        path vcf_tbi           // VCF 索引文件 (.tbi)
        path fasta             // 参考基因组 FASTA
        path fasta_fai         // 参考基因组索引 (.fai)
        val genome_assembly    // 基因组版本: GRCh37 或 GRCh38
        val upstream_distance  // 上游基因距离阈值 bp (默认 5000)
        val downstream_distance // 下游基因距离阈值 bp (默认 5000)
        val cache_dir          // VEP 缓存目录 (可选)
        val extra_args         // 额外参数 (可选)

    output:
        path "*.vep.mei.vcf.gz", emit: vep_vcf
        path "*.vep.mei.vcf.gz.tbi", emit: vep_vcf_tbi
        path "*.vep.mei.summary.txt", emit: summary

    script:
    def assembly = genome_assembly ?: 'GRCh38'
    def upstream = upstream_distance ?: 5000
    def downstream = downstream_distance ?: 5000
    def cache = cache_dir ? "--cache --dir_cache ${cache_dir}" : '--cache'
    def extra = extra_args ?: ''
    def sample_id = vcf.baseName.replaceAll(/\.(te\.result\.vcf|vcf|vcf\.gz)$/, '')
    """
    # 设置 VEP 缓存目录
    if [ "${cache_dir}" != "null" ] && [ -d "${cache_dir}" ]; then
        VEP_CACHE="${cache_dir}"
    else
        VEP_CACHE="/opt/vep/.vep"
    fi

    # 运行 VEP MEI 注释
    # --distance: 设置上下游基因距离阈值，MEI 影响范围较大
    # --symbol: 输出基因符号
    # --canonical: 标注 canonical 转录本
    # --domains: 输出蛋白结构域信息
    vep \\
        --input_file ${vcf} \\
        --output_file ${sample_id}.vep.mei.vcf \\
        --format vcf \\
        --vcf \\
        --assembly ${assembly} \\
        ${cache} \\
        --offline \\
        --dir_cache \${VEP_CACHE} \\
        --fasta ${fasta} \\
        --distance ${upstream},${downstream} \\
        --symbol \\
        --canonical \\
        --domains \\
        --protein \\
        --force_overwrite \\
        --no_stats \\
        ${extra}

    # 压缩 VCF 文件并创建索引
    bgzip -c ${sample_id}.vep.mei.vcf > ${sample_id}.vep.mei.vcf.gz
    tabix -p vcf ${sample_id}.vep.mei.vcf.gz

    # 生成 MEI 注释摘要
    echo "VEP MEI Annotation Summary for ${sample_id}" > ${sample_id}.vep.mei.summary.txt
    echo "Input VCF: ${vcf}" >> ${sample_id}.vep.mei.summary.txt
    echo "Output VCF: ${sample_id}.vep.mei.vcf.gz" >> ${sample_id}.vep.mei.summary.txt
    echo "Assembly: ${assembly}" >> ${sample_id}.vep.mei.summary.txt
    echo "Upstream/Downstream distance: ${upstream}/${downstream} bp" >> ${sample_id}.vep.mei.summary.txt
    echo "" >> ${sample_id}.vep.mei.summary.txt
    echo "Total MEI events annotated:" >> ${sample_id}.vep.mei.summary.txt
    zcat ${sample_id}.vep.mei.vcf.gz | grep -v "^#" | wc -l >> ${sample_id}.vep.mei.summary.txt
    echo "" >> ${sample_id}.vep.mei.summary.txt
    echo "MEI types breakdown:" >> ${sample_id}.vep.mei.summary.txt
    zcat ${sample_id}.vep.mei.vcf.gz | grep -v "^#" | grep -o "RN=[A-Za-z0-9_-]*" | sort | uniq -c >> ${sample_id}.vep.mei.summary.txt
    echo "" >> ${sample_id}.vep.mei.summary.txt
    echo "Consequence types:" >> ${sample_id}.vep.mei.summary.txt
    zcat ${sample_id}.vep.mei.vcf.gz | grep -v "^#" | cut -f8 | grep -o "Consequence=[A-Za-z_]*" | sed 's/Consequence=//' | sort | uniq -c >> ${sample_id}.vep.mei.summary.txt
    """
}