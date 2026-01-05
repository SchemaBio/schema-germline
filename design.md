# Schema Germline 系统设计文档

## 1. 项目背景

### 1.1 目标
构建一个高效、可扩展的外显子测序（WES）分析系统，支持：
- 单样本全流程分析
- 家系（Trio/多人）基因型叠加分析
- 前端友好的任务管理接口

### 1.2 设计原则
- **分离关注点**：流程引擎与分析逻辑分离
- **增量处理**：新样本入库后自动纳入家系分析
- **可扩展性**：支持本地/HPC/云端部署

---

## 2. 系统架构

```
┌─────────────────────────────────────────────────────────────────┐
│                         前端系统                                  │
│         (Web UI / 任务提交 / 结果展示 / 家系管理)                 │
└─────────────────────────────┬───────────────────────────────────┘
                              │ HTTP REST API
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      Python API 服务层                           │
│  ┌─────────────────────────────────────────────────────────────┐│
│  │ 任务管理 API                                                 ││
│  │   - 样本入库 / 状态查询 / 结果获取                           ││
│  │   - 任务提交 / 取消 / 重试                                   ││
│  └─────────────────────────────────────────────────────────────┘│
│  ┌─────────────────────────────────────────────────────────────┐│
│  │ 家系分析引擎                                                 ││
│  │   - gVCF 库管理                                              ││
│  │   - 家系成员检查                                             ││
│  │   - GLNexus 联合分型                                         ││
│  │   - Slivar 家系过滤                                          ││
│  └─────────────────────────────────────────────────────────────┘│
│  ┌─────────────────────────────────────────────────────────────┐│
│  │ 分析结果缓存                                                 ││
│  │   - 家系联合 VCF                                             ││
│  │   - 分析结果（TSV/JSON）                                     ││
│  └─────────────────────────────────────────────────────────────┘│
└─────────────────────────────┬───────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                   Nextflow 核心分析流程                          │
│  ┌─────────────────────────────────────────────────────────────┐│
│  │ wes_single.nf                                                ││
│  │   - FASTP 质控                                               ││
│  │   - BWA/BWA-MEM2 比对                                       ││
│  │   - GATK MarkDuplicates                                      ││
│  │   - DeepVariant 变异检测                                     ││
│  │   - WhatsHap 定相                                            ││
│  │   - VEP 注释                                                 ││
│  │   - GenMod 遗传模式                                          ││
│  │   - Slivar 过滤                                              ││
│  │   - CNV/STR/MT 专项分析                                      ││
│  └─────────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
              ┌───────────────┴───────────────┐
              │      HPC 调度器 / 本地         │
              │  (SLURM / LSF / K8s / local)  │
              └───────────────────────────────┘
```

---

## 3. 组件说明

### 3.1 Python API 服务层

| 模块 | 功能 |
|------|------|
| **任务管理** | 样本状态追踪、Nextflow 任务提交/监控 |
| **家系分析** | gVCF 合并、联合分型、家系特异性过滤 |
| **结果服务** | VCF/TSV 导出、结果缓存管理 |

#### API 接口设计

```
POST /api/samples/{sample_id}/analyze    # 提交样本分析
GET  /api/samples/{sample_id}/status     # 查询分析状态
GET  /api/samples/{sample_id}/result     # 获取分析结果
GET  /api/samples/{sample_id}/gvcf       # 获取 gVCF

POST /api/families/{family_id}/analyze   # 触发家系分析
GET  /api/families/{family_id}/result    # 获取家系分析结果
GET  /api/families/{family_id}/members   # 查询家系成员
```

### 3.2 Nextflow 分析流程

#### wes_single.nf 流程

```
输入：FASTQ 文件（read1, read2）
输出：
├── BAM/CRAM       # 比对结果
├── VCF            # 变异检测结果
├── g.VCF          # gVCF（用于家系分析）
├── 质控报告       # fastp, coverage, sex_check
└── 过滤结果       # Slivar 过滤后的变异
```

#### 关键输出

| 文件 | 用途 |
|------|------|
| `{sample_id}.g.vcf.gz` | 家系联合分型 |
| `{sample_id}.deepvariant.vcf.gz` | 原始变异 |
| `{sample_id}.slivar.vcf.gz` | 过滤后变异 |

### 3.3 gVCF 库设计

```
/data/
├── gvcf/                          # 单样本 gVCF 库
│   ├── {sample_id}/
│   │   ├── {sample_id}.g.vcf.gz
│   │   └── {sample_id}.g.vcf.gz.tbi
│   └── ...
│
├── families/                      # 家系联合 VCF（按需生成）
│   ├── {family_id}.joint.vcf.gz
│   └── {family_id}.joint.vcf.gz.tbi
│
└── cache/                         # 分析结果缓存
    ├── {family_id}_denovo.vcf.gz
    ├── {family_id}_AR.vcf.gz
    └── ...
```

---

## 4. 家系分析流程

### 4.1 样本入库

```python
def add_sample_to_library(sample_gvcf, sample_info):
    """
    1. gVCF 入库
    2. 检查家系完整性
    3. 如完整，触发家系预合并
    """
    # 保存 gVCF
    save_gvcf(sample_gvcf, sample_info.sample_id)

    # 检查家系
    family_members = get_family_members(sample_info.family_id)
    if all_members_complete(family_members):
        # 自动触发家系合并
        merge_family(sample_info.family_id)
```

### 4.2 家系联合分型

```python
def merge_family(family_id):
    """
    使用 GLNexus 合并家系成员 gVCF
    """
    gvcfs = get_family_gvcfs(family_id)

    # GLNexus 联合分型
    joint_vcf = run_glnexus(gvcfs)

    # 保存
    save_family_joint(family_id, joint_vcf)
```

### 4.3 家系特异性分析

```python
def analyze_family(family_id, proband_id, analysis_type):
    """
    使用 Slivar 进行家系特异性分析
    """
    joint_vcf = get_family_joint(family_id)
    ped_file = get_ped_file(family_id)

    # 家系过滤
    if analysis_type == 'denovo':
        expr = 'denovo'
    elif analysis_type == 'AR':
        expr = 'recessive'
    elif analysis_type == 'Xlinked':
        expr = 'x_linked'

    result = slivar_filter(joint_vcf, ped_file, expr)

    # 提取先证者结果
    return extract_sample(result, proband_id)
```

---

## 5. 部署方案

### 5.1 运行环境

| 组件 | 要求 |
|------|------|
| Python API | Python 3.10+, FastAPI/Flask |
| Nextflow | Java 11+, Nextflow 23.04+ |
| HPC | SLURM/LSF/K8s |
| 存储 | 共享文件系统（NFS/Lustre） |

### 5.2 容器化

```dockerfile
# API 服务
FROM python:3.10-slim
RUN pip install fastapi uvicorn

# Nextflow 流程（已有）
# ghcr.io/pzweuj/mapping:2025dec
```

### 5.3 资源规划

| 组件 | CPU | 内存 |
|------|-----|------|
| API 服务 | 4 | 8GB |
| Nextflow 任务 | 动态 | 动态 |
| 家系合并 (GLNexus) | 8 | 32GB |

---

## 6. 数据流

### 6.1 单样本分析

```
用户提交
    │
    ▼
Python API → 提交 Nextflow 任务 → 执行 wes_single.nf
                                        │
                                        ▼
                                  完成后回调 API
                                        │
                                        ▼
                                  gVCF 入库
                                        │
                                        ▼
                                  返回结果 URL
```

### 6.2 家系分析

```
用户选择家系 → 触发家系分析
                    │
                    ▼
            检查 gVCF 是否完整
                    │
        ┌───────────┴───────────┐
        ▼                       ▼
     不完整                  完整
     (等待)            GLNexus 合并
                            │
                            ▼
                    Slivar 过滤
                            │
                            ▼
                    返回结果
```

---

## 7. 扩展规划

### 7.1 短期
- [ ] Python API 服务实现
- [ ] Nextflow 流程容器化
- [ ] 基础家系分析功能

### 7.2 中期
- [ ] 结果缓存机制优化
- [ ] 增量 gVCF 合并
- [ ] 分布式任务调度

### 7.3 长期
- [ ] 多家系联合分析
- [ ] 群体频率计算
- [ ] 机器学习辅助筛选

---

## 8. 目录结构

```
schema-germline/
├── main.nf                    # Nextflow 入口
├── wes_single.nf              # 单样本流程
├── design.md                  # 本文档
├── conf/
│   └── modules.config         # 模块配置
├── bin/                       # 辅助脚本
├── modules/
│   └── local/                 # Nextflow 模块
└── workflows/
    └── cnv_baseline.nf       # CNV 基线构建
```

---

## 9. 参考资料

- Nextflow: https://www.nextflow.io/
- DeepVariant: https://github.com/google/deepvariant
- GLnexus: https://github.com/dnanexus/GLnexus
- Slivar: https://github.com/brentp/slivar
