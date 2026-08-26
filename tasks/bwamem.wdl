version 1.2

task BwaAlign {
    input {
        String prefix
        File read_1
        File read_2
        Directory ref_dir
        String ref_fasta_name
        Int threads
    }

    Int actual_threads = if threads > 32 then 32 else threads
    Int memory_gb = actual_threads * 2

    command <<<
        set -e
        
        # 默认使用标准 bwa
        USE_MINIBWA=false
        if [ ~{actual_threads} -ge 16 ]; then
            echo "[INFO] Threads >= 16, probing for minibwa acceleration..."
            if grep -q -E "avx2|avx512" /proc/cpuinfo; then
                echo "[INFO] CPU supports AVX2/AVX512."
                MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
                if [ "$MEM_KB" -gt 16777216 ]; then
                    echo "[INFO] Memory is sufficient ($MEM_KB KB)."
                    if [ -f "~{ref_dir}/~{ref_fasta_name}.l2b" ] && \
                       [ -f "~{ref_dir}/~{ref_fasta_name}.mbw" ]; then
                        echo "[INFO] minibwa indexes found!"
                        USE_MINIBWA=true
                    else
                        echo "[WARN] minibwa index missing (.l2b or .mbw). Falling back to standard bwa."
                    fi
                else
                    echo "[WARN] Insufficient memory for minibwa. Falling back."
                fi
            else
                echo "[WARN] CPU lacks AVX2/AVX512. Falling back."
            fi
        else
            echo "[INFO] Threads < 16. Using standard bwa."
        fi

        # 注意：无论是谁跑，最终吐出的 BAM 名字必须一样，保证下游无缝衔接
        
        if [ "$USE_MINIBWA" = true ]; then
            echo "[RUNNING] *** minibwa map ***"
            minibwa map \
                -R "@RG\tID:~{prefix}\tSM:~{prefix}\tPL:SchemaBio\tPU:Germline" \
                -t ~{actual_threads} \
                ~{ref_dir}/~{ref_fasta_name} \
                ~{read_1} ~{read_2} | \
            samtools sort -@ 4 -m 2G -o ~{prefix}.sorted.bam -

            # 原 bwa-mem2 命令暂时保留，便于后续恢复：
            # bwa-mem2 mem \
            #     -t ~{actual_threads} \
            #     -M -R "@RG\tID:~{prefix}\tSM:~{prefix}\tPL:SchemaBio\tPU:Germline" \
            #     ~{ref_dir}/~{ref_fasta_name} \
            #     ~{read_1} ~{read_2} | \
            # samtools sort -@ 4 -m 2G -o ~{prefix}.sorted.bam -
        else
            echo "[RUNNING] *** bwa mem ***"
            bwa mem \
                -t ~{actual_threads} \
                -M -R "@RG\tID:~{prefix}\tSM:~{prefix}\tPL:SchemaBio\tPU:Germline" \
                ~{ref_dir}/~{ref_fasta_name} \
                ~{read_1} ~{read_2} | \
            samtools sort -@ 4 -m 2G -o ~{prefix}.sorted.bam -
        fi
        
        # 统一建索引
        samtools index ~{prefix}.sorted.bam
    >>>

    output {
        File out_bam = "~{prefix}.sorted.bam"
        File out_bai = "~{prefix}.sorted.bam.bai"
    }

    runtime {
        cpu: actual_threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/mapping:v1.1.0"
    }
}
