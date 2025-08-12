from Bio import SeqIO

input_file = "gencode.v44.transcripts.fa"
output_file = "transcripts.fasta"
num_records = 1000

records = SeqIO.parse(input_file, "fasta")
subset = (record for i, record in enumerate(records) if i < num_records)
SeqIO.write(subset, output_file, "fasta")

print(f"前 {num_records} 筆 transcript 已儲存至 {output_file}")