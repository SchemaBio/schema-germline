### 下载hg38参考基因组

使用ENSEMBL提供的版本

https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

注意使用primary_assembly而不是top_level版本，primary_assembly不包含多态性Contig，在分析上效果更好。

#### 建立BWA索引

大内存机器建议使用BWA-MEM2，使用BWA-MEM2建立索引

```
bwa-mem2 index Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

小内存机器使用BWA建立索引

```
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

#### 建立samtools索引

```
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```