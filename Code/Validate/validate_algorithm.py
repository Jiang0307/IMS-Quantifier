import os
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import itertools
import gzip

# ---------- 參數 ----------
test_num = 10000000
k = 7
short_read_length = 100
alphabet = ['A', 'C', 'G', 'T']

# ---------- 資料夾路徑 ----------
transcript_folder = r"D:\Data\Transcript"
read_folder = r"D:\Data\Read"

# ---------- 讀取 fasta transcripts ----------
def load_transcripts(folder_path):
    transcripts = []
    for file in os.listdir(folder_path):
        if file.endswith(".fasta") or file.endswith(".fa"):
            file_path = os.path.join(folder_path, file)
            for record in SeqIO.parse(file_path, "fasta"):
                transcripts.append(str(record.seq).upper())
    return transcripts

# ---------- 讀取 fastq reads（支援 .gz）----------
def load_reads(folder_path):
    reads = []
    for file in os.listdir(folder_path):
        if file.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            file_path = os.path.join(folder_path, file)
            open_func = gzip.open if file.endswith(".gz") else open
            with open_func(file_path, "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    seq = str(record.seq).upper()
                    if len(seq) >= k:
                        reads.append(seq)
    return reads

# ---------- Bitvector 建立函數 ----------
def build_bitvector(sequence, k, kmer_to_index, num_kmers):
    bitvector = np.zeros(num_kmers, dtype=bool)

    # ✅ 若有非 ACGT 字元，標記為 discard
    if any(base not in {'A', 'C', 'G', 'T'} for base in sequence):
        discard_message_buffer.append(f"Current read contains non-ATCG bp, will discard result. Read sequence : {sequence}")
        return bitvector, True  # True 表示應跳過該 read

    for i in range(len(sequence) - k + 1):
        k_mer = sequence[i:i+k]
        idx = kmer_to_index[k_mer]
        bitvector[idx] = True

    return bitvector, False

# ---------- 讀取資料 ----------
transcripts = load_transcripts(transcript_folder)
reads = load_reads(read_folder)

num_transcripts = len(transcripts)
num_reads = len(reads)

if num_transcripts == 0 or num_reads == 0:
    print("Cannot find transcripts or reads, please check the folder")
    exit()

# ---------- 建所有可能的 k-mers ----------
all_kmers = [''.join(p) for p in itertools.product(alphabet, repeat=k)]
kmer_to_index = {kmer: idx for idx, kmer in enumerate(sorted(all_kmers))}
num_kmers = len(all_kmers)

# ---------- 建 index table ----------
index_table = defaultdict(lambda: np.zeros(num_transcripts, dtype=bool))
for idx, transcript in enumerate(transcripts):
    for i in range(len(transcript) - k + 1):
        k_mer = transcript[i:i+k]
        index_table[k_mer][idx] = True

# ---------- transcript bitvectors ----------
transcript_bitvectors = \
[
    build_bitvector(transcript, k, kmer_to_index, num_kmers)[0]
    for transcript in transcripts
]

# ---------- 比對每條 read ----------
inconsistent_cases = []
discard_reads = []
discard_message_buffer = []  # 暫存被 discard 的訊息
counter = 0

for read_id, short_read in enumerate(reads):
    if counter < test_num:
        if len(short_read) < k:
            continue

        # ---------- 如果 read 含非ATCG，跳過 ----------
        read_bv, contains_non_acgt = build_bitvector(short_read, k, kmer_to_index, num_kmers)
        if contains_non_acgt:
            discard_reads.append(read_id)
            continue

        # ---------- kallisto 比對 ----------
        result = np.ones(num_transcripts, dtype=bool)
        for i in range(len(short_read) - k + 1):
            k_mer = short_read[i:i+k]
            if k_mer in index_table:
                result = np.logical_and(result, index_table[k_mer])
            else:
                result = np.zeros(num_transcripts, dtype=bool)
                break
        kallisto_set = set(np.where(result == True)[0])

        # ---------- Bitvector 比對 ----------
        read_positions = np.where(read_bv == True)[0]
        bv_set = set()
        for t_idx, t_bv in enumerate(transcript_bitvectors):
            if all(t_bv[i] for i in read_positions):
                bv_set.add(t_idx)

        # ---------- 檢查一致性 ----------
        if kallisto_set != bv_set:
            inconsistent_cases.append(read_id)
            print(f"\n[DISTINCT] Read ID: {read_id}")
            print(f" - kallisto equivalence class  : {sorted(kallisto_set)}")
            print(f" - bitvector equivalence class : {sorted(bv_set)}")

        counter += 1

# ---------- 統計輸出 ----------
print(f"\n - k = {k}")
print(f" - Total {(counter)} reads tested (excluding discarded)")
print(f" - Identical result count ： {(counter - len(inconsistent_cases))}")
print(f" - Distinct result count : {len(inconsistent_cases)}")
print(f" - Distinct read ID ： {inconsistent_cases}")
print(f" - Discarded read count：{len(discard_reads)}")

# if discard_message_buffer:
#     print("\n[Discarded Reads Summary]")
#     for line in discard_message_buffer:
#         print(line)