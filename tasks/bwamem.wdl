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
    # Reserve enough memory for minibwa even when using fewer alignment threads.
    Int memory_gb = if actual_threads < 16 then 32 else actual_threads * 2

    command <<<
        set -euo pipefail
        
        # minibwa 为首选；线程数不决定比对工具，索引、CPU 和内存决定可用性。
        USE_MINIBWA=false
        if command -v minibwa >/dev/null 2>&1 && \
           [ -r "~{ref_dir}/~{ref_fasta_name}.l2b" ] && [ -s "~{ref_dir}/~{ref_fasta_name}.l2b" ] && \
           [ -r "~{ref_dir}/~{ref_fasta_name}.mbw" ] && [ -s "~{ref_dir}/~{ref_fasta_name}.mbw" ] && \
           grep -q -E "avx2|avx512" /proc/cpuinfo; then
            MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
            if [ "$MEM_KB" -gt 16777216 ]; then
                USE_MINIBWA=true
                echo "[INFO] Using preferred minibwa indexes."
            else
                echo "[WARN] Insufficient memory for minibwa ($MEM_KB KB)."
            fi
        fi
        if [ "$USE_MINIBWA" = false ]; then
            echo "[WARN] minibwa is unavailable on this node; checking standard BWA fallback."
            command -v bwa >/dev/null 2>&1 || { echo "[ERROR] Neither minibwa nor standard BWA is usable." >&2; exit 1; }
            for suffix in amb ann bwt pac sa; do
                [ -r "~{ref_dir}/~{ref_fasta_name}.$suffix" ] && [ -s "~{ref_dir}/~{ref_fasta_name}.$suffix" ] || \
                    { echo "[ERROR] minibwa is unavailable and standard BWA index .$suffix is missing or empty." >&2; exit 1; }
            done
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
