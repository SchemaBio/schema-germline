# Schema Germline Pipeline

基于 WDL 的全外显子胚系变异检测流程。

## 快速开始

### 1. 环境要求

- Docker
- miniwdl (`pip install miniwdl`)
- 自部署镜像已拉取（见下文）

#### HPC
- Slurm (`pip install miniwdl-slurm`)
- LSF (`pip install miniwdl-lsf`)

### 2. 拉取镜像

```bash
docker pull docker.schema-bio.com/schemabio/gatk:4.6.2.0
docker pull docker.schema-bio.com/schemabio/deepvariant:1.10.0
docker pull docker.schema-bio.com/schemabio/vep:115.2
docker pull docker.schema-bio.com/schemabio/whatshap:2.8
docker pull docker.schema-bio.com/schemabio/mapping:v1.0.0
```

### 3. 准备配置文件

每个样本一个 JSON 文件，如 `inputs/single.json`:

```json
{
    "SingleWES.prefix": "proband",
    "SingleWES.read_1": "/mnt/data/test/sample_R1.fq.gz",
    "SingleWES.read_2": "/mnt/data/test/sample_R2.fq.gz",
    "SingleWES.fasta": "/database/hg38.fa",
    "SingleWES.bed": "/mnt/data/test/capture.bed",
    "SingleWES.flank_size": 50,
    "SingleWES.assembly": "GRCh38",
    "SingleWES.ref_dir": "/database",
    "SingleWES.cache_dir": "/database/vep",
    "SingleWES.schema_bundle": "/database/schema_bundle",
    "SingleWES.tmp_dir": "/tmp_workspace"
}
```

> bwa 索引文件需与 fasta 同目录同前缀，无需单独配置。

### 4. 配置 miniwdl

在 `conf/local.cfg` 配置 Docker 挂载和资源。HPC环境请查看conf下对应的配置文件。

### 5. 运行流程

```bash
miniwdl run single.wdl \
    --cfg conf/local.cfg \
    -i inputs/single.json \
    --dir /mnt/data/output \
    -p /path/to/schema-germline
```

## License

MIT