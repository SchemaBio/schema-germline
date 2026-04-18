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
//   - custom: 自定义数据库注释 (gnomAD, AlphaMissense, etc.)
//   - plugin: VEP 插件 (FlankingSequence, AnnotateClinVar, MissenseZscoreTranscript)
//
// 支持的自定义数据库 (通过 --custom 参数添加):
//   - gnomAD: 人群频率数据库 (AF, AC, AN 等字段)
//   - AlphaMissense: 变异致病性预测评分
//   - EVOScore2: 进化保守性评分
//   - Pangolin: 基因组注释
//
// 支持的插件 (通过 --plugin 参数添加):
//   - FlankingSequence: 输出变异位点上下游序列
//   - AnnotateClinVar: ClinVar 临床意义注释 (CLNSIG, CLNDN, CLNSTAR)
//   - MissenseZscoreTranscript: 转录本级别 Missense Z-score

process VEP_ANNOTATE {
    tag "VEP on ${vcf.baseName}"
    label 'process_medium'
    label 'vep'

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
        // 自定义数据库文件 (可选)
        path gnomad_vcf        // gnomAD 人群频率数据库 VCF
        path gnomad_tbi        // gnomAD 索引文件
        path alphamissense_vcf // AlphaMissense 变异致病性预测 VCF
        path alphamissense_tbi // AlphaMissense 索引文件
        path evoscore_vcf      // EVOScore2 进化保守性评分 VCF
        path evoscore_tbi      // EVOScore 索引文件
        path pangolin_vcf      // Pangolin 基因组注释 VCF
        path pangolin_tbi      // Pangolin 索引文件
        // 插件相关参数 (可选)
        val flanking_seq_len   // FlankingSequence 插件: 上下游序列长度 (默认 10)
        path clinvar_vcf       // AnnotateClinVar 插件: ClinVar VCF 文件
        path clinvar_tbi       // ClinVar 索引文件
        path missense_bed      // MissenseZscoreTranscript 插件: Missense Z-score BED 文件

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

    // 构建自定义数据库参数
    def custom_args = ''
    if (gnomad_vcf) {
        custom_args += " --custom ${gnomad_vcf},gnomad,vcf,exact,force_overwrite"
    }
    if (alphamissense_vcf) {
        custom_args += " --custom ${alphamissense_vcf},AlphaMissense,vcf,exact,force_overwrite"
    }
    if (evoscore_vcf) {
        custom_args += " --custom ${evoscore_vcf},EVOScore2,vcf,exact,force_overwrite"
    }
    if (pangolin_vcf) {
        custom_args += " --custom ${pangolin_vcf},Pangolin,vcf,exact,force_overwrite"
    }

    // 构建插件参数
    def plugin_args = ''
    if (flanking_seq_len) {
        plugin_args += " --plugin FlankingSequence,${flanking_seq_len}"
    }
    if (clinvar_vcf) {
        plugin_args += " --plugin AnnotateClinVar,clinvar_file=${clinvar_vcf},fields=CLNSIG,CLNDN,CLNSTAR"
    }
    if (missense_bed) {
        plugin_args += " --plugin MissenseZscoreTranscript,${missense_bed},1"
    }
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
    # --custom 参数添加自定义数据库注释:
    #   格式: --custom <file>,<name>,<format>,<type>,<force_overwrite>
    #   gnomad: 人群频率 (AF, AC, AN, Homo, etc.)
    #   AlphaMissense: 致病性预测评分 (am_pathogenicity, am_class)
    #   EVOScore2: 进化保守性评分
    #   Pangolin: 基因组注释
    # --plugin 参数添加插件注释:
    #   FlankingSequence: 输出变异位点上下游序列
    #   AnnotateClinVar: ClinVar 临床意义 (CLNSIG, CLNDN, CLNSTAR)
    #   MissenseZscoreTranscript: 转录本级别 Missense Z-score
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
        ${custom_args} \\
        ${plugin_args} \\
        ${extra}

    # 压缩 VCF 文件并创建索引
    bgzip -c ${sample_id}.vep.vcf > ${sample_id}.vep.vcf.gz
    tabix -p vcf ${sample_id}.vep.vcf.gz

    # 生成注释摘要
    echo "VEP Annotation Summary for ${sample_id}" > ${sample_id}.vep.summary.txt
    echo "Input VCF: ${vcf}" >> ${sample_id}.vep.summary.txt
    echo "Output VCF: ${sample_id}.vep.vcf.gz" >> ${sample_id}.vep.summary.txt
    echo "Assembly: ${assembly}" >> ${sample_id}.vep.summary.txt
    echo "" >> ${sample_id}.vep.summary.txt
    echo "Custom databases used:" >> ${sample_id}.vep.summary.txt
    if [ -n "${gnomad_vcf}" ]; then echo "  - gnomAD: ${gnomad_vcf}" >> ${sample_id}.vep.summary.txt; fi
    if [ -n "${alphamissense_vcf}" ]; then echo "  - AlphaMissense: ${alphamissense_vcf}" >> ${sample_id}.vep.summary.txt; fi
    if [ -n "${evoscore_vcf}" ]; then echo "  - EVOScore2: ${evoscore_vcf}" >> ${sample_id}.vep.summary.txt; fi
    if [ -n "${pangolin_vcf}" ]; then echo "  - Pangolin: ${pangolin_vcf}" >> ${sample_id}.vep.summary.txt; fi
    echo "" >> ${sample_id}.vep.summary.txt
    echo "Plugins used:" >> ${sample_id}.vep.summary.txt
    if [ -n "${flanking_seq_len}" ]; then echo "  - FlankingSequence: ${flanking_seq_len} bp" >> ${sample_id}.vep.summary.txt; fi
    if [ -n "${clinvar_vcf}" ]; then echo "  - AnnotateClinVar: ${clinvar_vcf}" >> ${sample_id}.vep.summary.txt; fi
    if [ -n "${missense_bed}" ]; then echo "  - MissenseZscoreTranscript: ${missense_bed}" >> ${sample_id}.vep.summary.txt; fi
    echo "" >> ${sample_id}.vep.summary.txt
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