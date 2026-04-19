# Schema Germline Pipeline

基于 Nextflow DSL2 的全外显子胚系变异检测流程。

## 快速开始

### 1. 环境要求

- Docker
- Nextflow
- 自部署镜像已拉取（见下文）

### 2. 拉取镜像

```bash
# 流程依赖镜像
docker pull docker.schema-bio.com/schemabio/germline:v1.0.0
docker pull docker.schema-bio.com/schemabio/mapping:v1.0.0
docker pull docker.schema-bio.com/schemabio/gatk:4.6.2.0
docker pull docker.schema-bio.com/schemabio/cnvkit:0.9.13
docker pull docker.schema-bio.com/schemabio/expansionhunter:5.0.0
docker pull docker.schema-bio.com/schemabio/deepvariant:1.10.0
docker pull docker.schema-bio.com/schemabio/vep:115.2
docker pull docker.schema-bio.com/schemabio/glnexus:v1.4.1
docker pull docker.schema-bio.com/schemabio/whatshap:2.8
docker pull docker.schema-bio.com/schemabio/tiea_wes:2.0.0
docker pull docker.schema-bio.com/schemabio/automap:1.3
docker pull docker.schema-bio.com/schemabio/peddy:0.4.8
```

### 3. 准备配置文件

每个样本一个 JSON 文件，如 `examples/single.json`:

```json
{
  "sample_id": "sample1",
  "read1": "/mnt/d/data/sample1_R1.fq.gz",
  "read2": "/mnt/d/data/sample1_R2.fq.gz",
  "reference": {
    "fasta": "/mnt/d/reference/hg38.fa"
  },
  "outdir": "/mnt/d/analysis/results"
}
```

> bwa/bwa-mem2 索引文件需与 fasta 同目录同前缀，无需单独配置。

### 4. 运行流程

```bash
nextflow run main.nf \
    --config examples/single.json \
    -profile local
```

## License

MIT