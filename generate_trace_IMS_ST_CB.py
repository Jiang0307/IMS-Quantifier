import math
import xml.etree.ElementTree as ET
from pathlib import Path
from bidict import bidict
from Code.utils import *

current_dir = Path(__file__).resolve().parent 
# ===================== pseudoalignment config ===================== #
k = 7
read_path = "D:\Data\Read"
transcript_path = "D:\Data\Transcript"
MAX_REQUESTS_PER_FILE = 150000

transcripts = load_transcripts(transcript_path)
all_reads = list(load_reads(read_path, k))
num_transcripts = len(transcripts)
num_reads = len(all_reads)
all_kmers_slim = extract_all_kmers(transcripts, all_reads, k)
kmer_to_index_slim = {kmer: idx for idx, kmer in enumerate(all_kmers_slim)}
num_kmers_slim = len(all_kmers_slim)
bitvector_length = len(all_kmers_slim)
print(f"[INFO] Bitvector length before compaction：{4**k}")
print(f"[INFO] Bitvector length after compaction ：{bitvector_length}")

# ===================== NAND flash config ===================== #
start_time = 0
read_request_interval = 220000
write_request_interval = 220000
write_counter = 0
read_counter = 0
write_request_counter = 0
read_request_counter = 0
MQSim_traces = []

# ===================== trace 儲存設定 ===================== #
ssdconfig_path = current_dir / "SSD_config.xml"
trace_dir = current_dir / "Traces"

if trace_dir.exists():
    for file in trace_dir.glob("MQSim_*.trace"):
        file.unlink()
else:
    trace_dir.mkdir(parents=True, exist_ok=True)

# ===================== SSD 架構設定讀取 ===================== #
tree = ET.parse(ssdconfig_path)
root = tree.getroot()
device_para_node = root.find("Device_Parameter_Set")
n_channels = int(device_para_node.find("Flash_Channel_Count").text)
n_chips_per_channel = int(device_para_node.find("Chip_No_Per_Channel").text)
flash_para_node = device_para_node.find("Flash_Parameter_Set")
n_dies_per_chip = int(flash_para_node.find("Die_No_Per_Chip").text)
n_planes_per_die = int(flash_para_node.find("Plane_No_Per_Die").text)
n_blocks_per_plane = int(flash_para_node.find("Block_No_Per_Plane").text)
n_pages_per_block = int(flash_para_node.find("Page_No_Per_Block").text)
page_capacity = int(flash_para_node.find("Page_Capacity").text)
sector_size = 512
n_sectors_per_page = int(page_capacity / sector_size)
word_lines = int(flash_para_node.find("Page_No_Per_Block").text)
reads_per_block = page_capacity * 8 # 每個block可容納的read數
PSU_independent_dies = n_channels * n_chips_per_channel * n_dies_per_chip

sectors_per_block = n_pages_per_block * n_sectors_per_page
sectors_per_die = n_blocks_per_plane * sectors_per_block
print(f"PSU_independent_dies : {PSU_independent_dies}")
# ===================== 計算寫入與讀取次數 ===================== #

# bitvector_length = 4 ** k
bits_per_block = word_lines / 2
blocks_per_bitvector = math.ceil(bitvector_length / bits_per_block)
read_groups = math.ceil(num_reads / reads_per_block)
total_block_required = read_groups * blocks_per_bitvector

if total_block_required >= PSU_independent_dies:
    print(f"total_block_required ≥ {PSU_independent_dies}")
    n_PSU = math.ceil(total_block_required / PSU_independent_dies)
    PSU_write_count = total_block_required
    PSU_read_count = math.ceil(num_transcripts * total_block_required)

else:
    print(f"total_block_required < {PSU_independent_dies}")
    n_PSU = 1
    PSU_write_count = total_block_required
    PSU_read_count = math.ceil(num_transcripts * total_block_required)

# ===================== WRITE =====================
for PSU in range(n_PSU):
    for channel in range(n_channels): # 這層for以下就是等於1 PSU
        for chip in range(n_chips_per_channel):
            for die in range(n_dies_per_chip):
                die_index = (channel * n_chips_per_channel * n_dies_per_chip) + (chip * n_dies_per_chip) + (die)
                block_start_lsa = (die_index * sectors_per_die) + (PSU * sectors_per_block)
                if write_counter < PSU_write_count: # 總共要寫PSU_write_count次
                    if write_counter % PSU_independent_dies == 0:
                        for i in range(n_pages_per_block): # 一次IMS寫入會寫滿一個block = 256個wordline(page)
                            page_start_lsa = block_start_lsa + i * (page_capacity/sector_size)
                            MQSim_trace = MQSim_Trace(start_time, 0, page_start_lsa, page_capacity/sector_size, 0) # Request_Size_In_Sectors = 1，NAND flash write時一樣會write整頁
                            MQSim_traces.append(MQSim_trace)
                            write_request_counter += 1
                            start_time += write_request_interval
                    write_counter += 1
                else:
                    break

# ===================== READ =====================
for PSU in range(n_PSU * num_transcripts):
    for channel in range(n_channels): # 這層for以下就是等於1 PSU
        for chip in range(n_chips_per_channel):
            for die in range(n_dies_per_chip):
                die_index = (channel * n_chips_per_channel * n_dies_per_chip) + (chip * n_dies_per_chip) + (die)
                block_start_lsa = (die_index * sectors_per_die) + (PSU * sectors_per_block)
                
                if read_counter < PSU_read_count: # 總共要搜尋PSU_read_count次
                    if read_counter % PSU_independent_dies == 0:
                        MQSim_trace = MQSim_Trace(start_time, 0, block_start_lsa, page_capacity/sector_size, 1) # Request_Size_In_Sectors = 1，NAND flash read時一樣會read整頁
                        MQSim_traces.append(MQSim_trace)
                        read_request_counter += 1
                    start_time += read_request_interval
                    read_counter += 1
                else:
                    break
    if read_counter >= PSU_read_count:
        break

# ===================== 輸出 MQSim_Trace =====================
for i in range(0, len(MQSim_traces), MAX_REQUESTS_PER_FILE):
    chunk = MQSim_traces[i:i + MAX_REQUESTS_PER_FILE]
    chunk_file = current_dir / "Traces" / f"MQSim_{i // MAX_REQUESTS_PER_FILE}.trace"
    output_trace(chunk_file, chunk)

num_trace_files = math.ceil(len(MQSim_traces) / MAX_REQUESTS_PER_FILE)
workload_output_path = current_dir / "Workload.xml"
generate_workload_xml(num_trace_files, workload_output_path)

print(f"\nnum_reads = {num_read_dict.inverse[num_reads]} , num_transcripts = {num_transcripts}")
print(f"k = {k} , bitvector_length = {bitvector_length} , blocks_per_bitvector = {blocks_per_bitvector}")
print(f"total_block_required = {total_block_required}")
print(f"write_requests = {math.ceil(write_request_counter)}")
print(f"read_requests = {math.ceil(read_request_counter)}")