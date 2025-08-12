import subprocess
import datetime
import re
import argparse
from pathlib import Path

# === 取得當前 .py 檔所在目錄 ===
current_dir = Path(__file__).resolve().parent

# === 解析參數 ===
parser = argparse.ArgumentParser()
parser.add_argument("-k", type=int, required=True, help="k-mer size (will load index_k.idx)")
args = parser.parse_args()

k_value = args.k

# === 組合 index 檔案完整路徑 ===
index_table_path = current_dir / "index_table" / f"index_{k_value}.idx"
output_path = current_dir / f"pseudoalignment_output_{k_value}"
fastq_path = Path("C:/Users/User/Desktop/GitHub/IMS-quantifier/Dataset/Read/SRR062634_1_subset.fastq")

# === 定義 kallisto 指令 ===
cmd = [
    str(current_dir / "kallisto.exe"),
    "pseudo",
    "-i", str(index_table_path),
    "-o", str(output_path),
    "--single",
    "-l", "100",
    "-s", "1",
    "-t", "8",
    str(fastq_path)
]

# === pattern for pseudoalignment 結束偵測 ===
pattern = r"\[quant\] processed [\d,]+ reads, [\d,]+ reads pseudoaligned"

print("[system] running:", " ".join(cmd))
start_time = datetime.datetime.now()
print(f"[system] start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

# === 執行 kallisto 並即時讀取輸出 ===
process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

pseudoaligned_logged = False

for line in process.stdout:
    print(line, end='')  # 即時印出 kallisto 輸出
    if not pseudoaligned_logged and re.search(pattern, line):
        end_time = datetime.datetime.now()
        duration = end_time - start_time
        duration_sec = duration.total_seconds()
        duration_us = int(duration_sec * 1000000)

        print(f"\n[system] end time: {end_time.strftime('%Y-%m-%d %H:%M:%S.%f')}")
        print(f"[system] duration: {duration_us} μs = {duration_sec:.6f} s")
        pseudoaligned_logged = True

process.wait()