version 1.2

task Whatshap {
    input {
        String prefix
        File bam
        File bai
        File vcf
        File vcf_tbi
        String fasta
        Int threads
        Directory ref_dir
    }

    Int memory_gb = threads * 2

    command <<<
        set -e # 遇到错误立刻停止

        # 1. 使用 tabix 提取 VCF 中包含变异的所有染色体/Contigs 名称
        tabix -l ~{vcf} > contigs.txt

        # 2. 创建临时目录存放分染色体的 VCF
        mkdir -p tmp_vcfs

        # 3. 使用 xargs 实现多线程并发处理各个染色体
        # -P ~{threads} 表示同时运行的进程数
        cat contigs.txt | xargs -I {} -P ~{threads} sh -c '
            whatshap phase \
                --chromosome {} \
                --reference=~{ref_dir}/~{fasta} \
                -o tmp_vcfs/~{prefix}.{}.phase.vcf \
                ~{vcf} \
                ~{bam}
        '

        # 4. 获取所有生成的子 VCF 列表并合并 (依赖 bcftools)
        # 将按染色体分散的 VCF 合并为一个完整的全基因组 VCF
        ls tmp_vcfs/~{prefix}.*.phase.vcf > vcf_list.txt
        bcftools concat -f vcf_list.txt -O z -o ~{prefix}.phase.vcf.gz
        
        # 5. 对最终结果建索引
        tabix -f -p vcf ~{prefix}.phase.vcf.gz
    >>>

    output {
        File out_vcf = "~{prefix}.phase.vcf.gz"
        File out_vcf_tbi = "~{prefix}.phase.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/whatshap:2.8"
    }
}