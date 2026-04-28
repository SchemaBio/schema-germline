version 1.2

task Peddy {
    input {
        String prefix
        File vcf
        File vcf_idx
        File ped
        String assembly
        Int threads
    }

    Int memory_gb = threads * 2

    command <<<
        set -e
    
        if [ "~{assembly}" == "GRCh38" ]; then
            fix_assembly="hg38"
        elif [ "~{assembly}" == "GRCh37" ]; then
            fix_assembly="hg19"
        else
            echo "Unsupported assembly: ~{assembly}" >&2
            exit 1
        fi
        
        # 将 VCF 和索引软链接到当前执行目录，确保前缀一致以便被顺利读取
        VCF_BASENAME=$(basename ~{vcf})
        ln -s ~{vcf} $VCF_BASENAME
        ln -s ~{vcf_idx} ${VCF_BASENAME}.tbi 
        # 如果你的索引是 .csi 格式，请将上一行的 .tbi 改为 .csi

        echo "开始运行 peddy..."
        peddy \
            -p ~{threads} \
            --sites ${fix_assembly} \
            --plot \
            --prefix ~{prefix} \
            $VCF_BASENAME \
            ~{ped}
    >>>

    output {
        File ped_check = "~{prefix}.ped_check.csv"
        File het_check = "~{prefix}.het_check.csv"
        File sex_check = "~{prefix}.sex_check.csv"
    }

    runtime {
        docker: "docker.schema-bio.com/schemabio/peddy:0.4.8"
        cpu: threads
        memory: "~{memory_gb}G"
    }
}
