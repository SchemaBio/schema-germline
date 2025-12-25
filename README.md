

```
schema-germline/
├── assets/                     # [资产] 静态文件，如 MultiQC 模板、Bed 文件、JSON Schema
│   └── multiqc_config.yaml
│
├── bin/                        # [脚本] 存放 Python/R 胶水代码 (自动进入 $PATH)
│   ├── vcf_to_parquet.py       # <--- 核心：将 VCF 转为 Parquet 的脚本
│   └── qc_to_json.py           # <--- 核心：将 QC 结果清洗为前端可读的 JSON
│
├── conf/                       # [配置] 基础设施与参数的控制中心
│   ├── base.config             # 基础资源配置 (定义默认 CPU/Memory)
│   ├── modules.config          # <--- 核心：定义每个 Process 用哪个 Docker 镜像、什么参数
│   ├── infrastructure.config   # 云端配置 (AWS/阿里云 Spot 实例策略)
│   └── test.config             # CI/CD 测试用的微型数据集配置
│
├── containers/                 # [镜像] Docker 构建文件 (Infrastructure as Code)
│   ├── base/                   # 基础镜像 (Samtools, BWA 等)
│   │   └── Dockerfile
│   ├── gatk/                   # GATK 专用镜像
│   │   └── Dockerfile
│   └── python/                 # 数据处理镜像 (Pandas, PyArrow)
│       └── Dockerfile
│
├── modules/                    # [模块] 原子化的任务单元 (Process)
│   ├── local/                  # 我们自己写的私有模块
│   │   ├── bwa_mem2.nf         # <--- 改良版：直接输出 CRAM
│   │   ├── cram_qc.nf          # CRAM 格式的质控
│   │   └── vcf_conversion.nf   # 调用 bin/vcf_to_parquet.py
│   └── nf-core/                # 从 nf-core 社区抄作业的通用模块 (如 FastQC)
│
├── workflows/                  # [流程] 将模块串联起来的业务逻辑
│   └── germline.nf             # 主流程逻辑
│
├── main.nf                     # [入口] 程序的启动入口
├── nextflow.config             # [总控] 全局配置文件
└── README.md                   # [文档] 包含架构图和运行指南
```