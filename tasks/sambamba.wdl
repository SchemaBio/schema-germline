version 1.2

task SambambaMarkdup {
    input {
        String prefix
        File bam
        File bai
        Int threads
    }

    Int memory_gb = threads * 2

    command <<<
        set -e

        # 确保临时目录存在，防止运行时报错
        mkdir -p ./tmp

        sambamba markdup \
            --nthreads ~{threads} \
            --overflow-list-size 1000000 \
            --hash-table-size 1000000 \
            --compression-level 1 \
            ~{bam} \
            ~{prefix}.markdup.bam \
            --tmpdir=./tmp \
            2> ~{prefix}.sambamba.log

        samtools index ~{prefix}.markdup.bam
        
        # 使用 flagstat 获取标准化的统计数据，比解析工具 log 更安全
        samtools flagstat ~{prefix}.markdup.bam > ~{prefix}.flagstat

        STATS_FILE="~{prefix}.sambamba.stats"

        # 安全地提取数据：使用 || true 防止 set -e 导致流程崩溃
        TOTAL_READS=$(grep "in total" ~{prefix}.flagstat | awk '{print $1}' || true)
        DUPLICATES=$(grep "duplicates" ~{prefix}.flagstat | awk '{print $1}' || true)

        # 兜底机制：如果为空则赋值为 0，防止 Awk 报错
        TOTAL_READS=${TOTAL_READS:-0}
        DUPLICATES=${DUPLICATES:-0}

        # 使用 Awk 计算文库复杂度 (Lander-Waterman equation)
        awk -v N="$TOTAL_READS" -v D="$DUPLICATES" '
        BEGIN {
            U = N - D;

            # 边界情况处理
            if (N == 0 || U == 0) {
                C = 0;
            } else if (U >= N) {
                C = N;
            } else {
                # 二分查找逼近求值
                lower = U;
                upper = N * 1000.0;

                for (i = 0; i < 100; i++) {
                    C = (lower + upper) / 2.0;
                    U_est = C * (1.0 - exp(-N / C));

                    diff = U_est - U;
                    if (diff < 0) diff = -diff;

                    if (diff < 0.1) {
                        break;
                    } else if (U_est < U) {
                        lower = C;
                    } else {
                        upper = C;
                    }
                }
            }

            # 打印与 GATK Picard 兼容的报告
            printf "Total Reads:          %d\n", N
            printf "Duplicate Reads:      %d\n", D
            printf "Unique Reads:         %d\n", U
            printf "PERCENT_DUPLICATION:  %.4f\n", (N > 0 ? (D / N) : 0)
            printf "--------------------------------------\n"
            printf "ESTIMATED_LIBRARY_SIZE: %d\n", C
        }' > "$STATS_FILE"

        cat "$STATS_FILE"
    >>>

    output {
        File markdup_bam = "~{prefix}.markdup.bam"
        File markdup_bai = "~{prefix}.markdup.bam.bai"
        File sambamba_stats = "~{prefix}.sambamba.stats"
        File flagstat_log = "~{prefix}.flagstat" 
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "docker.schema-bio.com/schemabio/mapping:v1.0.0"
    }
}