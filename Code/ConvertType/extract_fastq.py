from Bio import SeqIO

# 讀取原始的FASTQ檔案
input_file = r"D:\Data\Read\ERR251006_10M_1.fastq"
output_file = "ERR251006_3200K.fastq"
num_reads_to_extract = 3200000  # 要取出的read數量

# 開啟輸出檔案
with open(output_file, "w") as out_f:
    # 記錄目前已取出的數量
    count = 0
    
    # 開啟原始的FASTQ檔案
    for record in SeqIO.parse(input_file, "fastq"):
        if count < num_reads_to_extract:
            SeqIO.write(record, out_f, "fastq")
            count += 1
        else:
            break