## Lamarck &nbsp; &nbsp; &nbsp; 2024-10-28
---
[RABV参考基因组的单基因链接](https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000859625.1/)
基因组序列数据集：`rabies_genomes.fasta`
模板基因序列数据集：`reference_n_gene.fasta`

---

**分析步骤**
1. 下载单基因的核苷酸序列文件：reference_n_gene.fasta；
2. 整理全基因组序列集：rabies_genomes.fasta；
3. 用全基因组序列集建库  makeblastdb；
4. 用单基因序列比对建好的库，输出为out文件；
5. 检查out文件中的Gene ID数目是否正确；
6. 根据out文件和全基因组序列集，提取单基因（1.py）；
7. 用mafft对单基因序列集进行多序列比对；
8. 用Aliview对多序列比对结果进行裁剪；

*建库*
```bash
makeblastdb -in rabies_genomes.fasta -dbtype nucl -out rabies_genome_db
```

*blast比对*
```bash
blastn -query reference_L_gene.fasta -db rabies_genome_db -out L_gene_results.out -outfmt 6 -evalue 1e-5 -max_target_seqs 2000 -word_size 7
```

*检查out文件中唯一Gene ID数目*
```bash
cut -f2 L_gene_results.out | sort | uniq | wc -l
```

*1.py 提取单基因序列*
```python
from Bio import SeqIO

# 读取比对结果
blast_results = "L_gene_results.out"
genome_file = "rabies_genomes.fasta"
output_file = "extracted_L_genes.fasta"

# 读取所有基因组序列
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

# 解析BLAST结果文件，提取G基因片段
extracted_sequences = []
with open(blast_results, "r") as f:
    for line in f:
        cols = line.strip().split("\t")
        sseqid = cols[1]
        sstart = int(cols[8])
        send = int(cols[9])

        # 获取并提取对应的G基因序列
        if sseqid in genome_dict:
            if sstart < send:
                g_gene_seq = genome_dict[sseqid].seq[sstart-1:send]
            else:
                g_gene_seq = genome_dict[sseqid].seq[send-1:sstart].reverse_complement()

            extracted_sequences.append(SeqIO.SeqRecord(g_gene_seq, id=sseqid, description=""))

# 将提取的G基因序列保存到输出文件
SeqIO.write(extracted_sequences, output_file, "fasta")
print(f"提取的L基因序列已保存到 {output_file}")

```
