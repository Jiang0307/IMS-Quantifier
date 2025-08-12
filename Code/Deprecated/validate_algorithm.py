import random
from collections import defaultdict
import numpy as np
import itertools

# ---------- 參數 ----------
num_transcripts = 1000
transcript_length = 1500
k = 5
short_read_length = 100
num_reads = 2000
alphabet = ['A', 'C', 'G', 'T']
insert_shared_prob = 0.05  # 每個位置插入 shared k-mer 的機率
num_shared_kmers = 200     # 共用 k-mer 數量

# ---------- 準備所有可能的 k-mers ----------
all_kmers = [''.join(p) for p in itertools.product(alphabet, repeat=k)]
kmer_to_index = {kmer: idx for idx, kmer in enumerate(sorted(all_kmers))}
num_kmers = len(all_kmers)

# ---------- 建立 shared k-mer pool ----------
shared_kmers_pool = random.sample(all_kmers, num_shared_kmers)

# ---------- 建構 transcript（帶有部分 shared k-mers） ----------
def generate_transcript(length, shared_pool, insert_prob):
    seq = []
    i = 0
    while i < length - k + 1:
        if random.random() < insert_prob:
            kmer = random.choice(shared_pool)
        else:
            kmer = ''.join(random.choices(alphabet, k=k))
        seq.append(kmer)
        i += 1
    # 組合回 DNA 序列
    return seq[0] + ''.join(kmer[-1] for kmer in seq[1:])

transcripts = [
    generate_transcript(transcript_length, shared_kmers_pool, insert_shared_prob)
    for _ in range(num_transcripts)
]

# ---------- 建 index table ----------
index_table = defaultdict(lambda: np.zeros(num_transcripts, dtype=bool))
for idx, transcript in enumerate(transcripts):
    for i in range(transcript_length - k + 1):
        k_mer = transcript[i:i+k]
        index_table[k_mer][idx] = True

# ---------- Bitvector 建立函數 ----------
def build_bitvector(sequence, k, kmer_to_index, num_kmers):
    bitvector = np.zeros(num_kmers, dtype=bool)
    for i in range(len(sequence) - k + 1):
        k_mer = sequence[i:i+k]
        idx = kmer_to_index.get(k_mer)
        if idx is not None:
            bitvector[idx] = True
    return bitvector

# ---------- transcript bitvectors ----------
transcript_bitvectors = [
    build_bitvector(transcript, k, kmer_to_index, num_kmers)
    for transcript in transcripts
]

# ---------- 測試多條 short reads ----------
inconsistent_cases = []

for read_id in range(num_reads):
    tid = random.randint(0, num_transcripts - 1)
    pos = random.randint(0, transcript_length - short_read_length)
    short_read = transcripts[tid][pos : pos + short_read_length]

    # ---------- Pseudo-code 比對 ----------
    result = np.ones(num_transcripts, dtype=bool)
    for i in range(len(short_read) - k + 1):
        k_mer = short_read[i:i+k]
        if k_mer in index_table:
            result = np.logical_and(result, index_table[k_mer])
        else:
            result = np.zeros(num_transcripts, dtype=bool)
            break

    # result[0] = 1
    # result[1] = 1
    # result[2] = 1
    # result[3] = 1
    # result[4] = 1

    pseudo_set = set(np.where(result == True)[0])

    # ---------- Bitvector 比對 ----------
    read_bv = build_bitvector(short_read, k, kmer_to_index, num_kmers)
    # read_bv[2] = 1
    # read_bv[3] = 1
    # read_bv[5] = 1
    # read_bv[7] = 1
    # read_bv[8] = 1
    
    read_positions = np.where(read_bv == True)[0]
    bv_set = set()
    for t_idx, t_bv in enumerate(transcript_bitvectors):
        if all(t_bv[i] for i in read_positions):
            bv_set.add(t_idx)

    # ---------- 比對結果是否一致 ----------
    if pseudo_set != bv_set:
        inconsistent_cases.append({
            "read_id": read_id,
            "source_tid": tid,
            "pseudo": sorted(pseudo_set),
            "bitvector": sorted(bv_set)
        })

    # ---------- 輸出每條 read 的結果 ----------
    # print(f"\nRead {read_id} | 來源 transcript: {tid}")
    # print(f"Pseudo-code 相容 transcripts : {sorted(pseudo_set)}")
    # print(f"Bitvector   相容 transcripts : {sorted(bv_set)}")

# ---------- 統計 ----------
print(f"\總共測試 {num_reads} 條 short read")
print(f"完全一致的比對數量：{num_reads - len(inconsistent_cases)}")
print(f"不一致的比對數量：{len(inconsistent_cases)}")