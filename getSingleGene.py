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
