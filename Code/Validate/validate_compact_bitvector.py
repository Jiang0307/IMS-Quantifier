import os
import itertools
import gzip
import numpy as np
from collections import defaultdict, Counter
from Bio import SeqIO

# ---------- 參數 ----------
test_num = 100
k = 13
alphabet = ['A', 'C', 'G', 'T']

# ---------- 資料夾路徑 ----------
transcript_folder = "D:\Data\Transcript"
read_folder = "D:\Data\Read"

# ---------- 讀取 transcripts ----------
def load_transcripts(folder_path):
    transcripts = []
    for file in os.listdir(folder_path):
        if file.endswith(".fasta") or file.endswith(".fa"):
            file_path = os.path.join(folder_path, file)
            for record in SeqIO.parse(file_path, "fasta"):
                transcripts.append(str(record.seq).upper())
    return transcripts

# ---------- 讀取 reads（用 generator + 限制筆數）----------
def load_reads(folder_path, k, max_reads=None):
    count = 0
    for file in os.listdir(folder_path):
        if file.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            file_path = os.path.join(folder_path, file)
            open_func = gzip.open if file.endswith(".gz") else open
            with open_func(file_path, "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    seq = str(record.seq).upper()
                    if len(seq) >= k:
                        yield seq
                        count += 1
                        if max_reads and count >= max_reads:
                            return

# ---------- 建 bitvector ----------
def build_bitvector(sequence, k, kmer_to_index, num_kmers):
    bitvector = np.zeros(num_kmers, dtype=bool)
    if any(base not in {'A', 'C', 'G', 'T'} for base in sequence):
        return bitvector, True
    for i in range(len(sequence) - k + 1):
        k_mer = sequence[i:i+k]
        idx = kmer_to_index.get(k_mer)
        if idx is not None:
            bitvector[idx] = True
    return bitvector, False

# ---------- 萃取實際 k-mers ----------
def extract_all_kmers(transcripts, reads, k):
    kmer_counter = Counter()
    for seq in transcripts + reads:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if set(kmer).issubset({'A', 'C', 'G', 'T'}):
                kmer_counter[kmer] += 1
    return sorted(kmer_counter.keys())

# ---------- 主程式 ----------
print("[STEP] Load transcripts")
transcripts = load_transcripts(transcript_folder)

print("[STEP] Load reads")
all_reads = list(load_reads(read_folder, k))

# ---------- 原始 k-mer index（full 4^k）----------
all_kmers_full = [''.join(p) for p in itertools.product(alphabet, repeat=k)]
kmer_to_index_full = {kmer: idx for idx, kmer in enumerate(sorted(all_kmers_full))}
num_kmers_full = len(all_kmers_full)
print(f"[INFO] Bitvector length before compaction：{num_kmers_full}")

# ---------- 壓縮版 k-mer index（只針對出現過的） ----------
# print("[STEP] Extracting k-mers")
all_kmers_slim = extract_all_kmers(transcripts, all_reads, k)
kmer_to_index_slim = {kmer: idx for idx, kmer in enumerate(all_kmers_slim)}
num_kmers_slim = len(all_kmers_slim)
print(f"[INFO] Bitvector length after compaction：{num_kmers_slim}")

# ---------- 建 index table（kallisto 模擬）----------
index_table = defaultdict(lambda: np.zeros(len(transcripts), dtype=bool))
for idx, transcript in enumerate(transcripts):
    for i in range(len(transcript) - k + 1):
        k_mer = transcript[i:i+k]
        index_table[k_mer][idx] = True

# ---------- transcript bitvector ----------
# transcript_bv_full = [
#     build_bitvector(t, k, kmer_to_index_full, num_kmers_full)[0]
#     for t in transcripts
# ]
transcript_bv_slim = [
    build_bitvector(t, k, kmer_to_index_slim, num_kmers_slim)[0]
    for t in transcripts
]

# ---------- 比對每條 read ----------
print(f"[STEP] 開始比對前 {test_num} 筆 reads")
inconsistent_cases = []
discard_reads = []
counter = 0

for read_id, short_read in enumerate(all_reads):
    if counter >= test_num:
        break
    if len(short_read) < k:
        continue

    # read_bv_full, _ = build_bitvector(short_read, k, kmer_to_index_full, num_kmers_full)
    read_bv_slim, contains_non_acgt = build_bitvector(short_read, k, kmer_to_index_slim, num_kmers_slim)
    if contains_non_acgt:
        discard_reads.append(read_id)
        continue

    # ---------- kallisto 比對 ----------
    result = np.ones(len(transcripts), dtype=bool)
    for i in range(len(short_read) - k + 1):
        k_mer = short_read[i:i+k]
        if k_mer in index_table:
            result = np.logical_and(result, index_table[k_mer])
        else:
            result = np.zeros(len(transcripts), dtype=bool)
            break
    kallisto_set = set(np.where(result == True)[0])

    # ---------- 原始 Bitvector 比對（已註解）----------
    # read_positions_full = np.where(read_bv_full == True)[0]
    # bv_full_set = set()
    # for t_idx, t_bv in enumerate(transcript_bv_full):
    #     if all(t_bv[i] for i in read_positions_full):
    #         bv_full_set.add(t_idx)

    # ---------- 壓縮 Bitvector 比對 ----------
    read_positions_slim = np.where(read_bv_slim == True)[0]
    bv_slim_set = set()
    for t_idx, t_bv in enumerate(transcript_bv_slim):
        if all(t_bv[i] for i in read_positions_slim):
            bv_slim_set.add(t_idx)

    # ---------- 結果比對 ----------
    if kallisto_set != bv_slim_set:
        inconsistent_cases.append(read_id)
        print(f"\n[DISTINCT] Read ID: {read_id}")
        print(f" - kallisto equivalence class  : {sorted(kallisto_set)}")
        # print(f" - bitvector (full) class      : {sorted(bv_full_set)}")
        print(f" - bitvector (slim) class      : {sorted(bv_slim_set)}")

    counter += 1

# ---------- 統計輸出 ----------
print(f"\n - k = {k}")
print(f" - Total reads tested         : {counter}")
print(f" - Identical result count     : {counter - len(inconsistent_cases)}")
print(f" - Distinct result count      : {len(inconsistent_cases)}")
print(f" - Distinct read ID list      : {inconsistent_cases}")
print(f" - Discarded read count       : {len(discard_reads)}")