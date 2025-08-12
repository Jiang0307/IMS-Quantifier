from Bio import SeqIO

input_fastq_path = r"D:\Data\Read\ERR251006_10M_1.fastq"
output_fa_path = r"D:\Data\Read\read_10M.fa"

with open(input_fastq_path, 'r') as fastq_file, open(output_fa_path, 'w') as fa_file:
    lines = fastq_file.readlines()
    seq_id = 0  # 初始化序列 ID 從 0 開始
    for i in range(0, len(lines), 4):
        # 提取 FASTQ 檔案中的序列
        sequence = lines[i+1].strip()  # 獲取序列
        fa_file.write(f">{seq_id}\n")  # 寫入 ID，ID 以 '>' 開頭
        fa_file.write(f"{sequence}\n")  # 寫入序列
        seq_id += 1  # ID 遞增