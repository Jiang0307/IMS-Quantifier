import math
import xml.etree.ElementTree as ET
from pathlib import Path

# ===================== 基本參數設定 ===================== #
perellelism_level = "channel"
k = 5
num_reads = 1000000
num_transcripts = 1000
word_lines = 256
reads_per_block = 131072  # 每個 block 可容納的 read 數
blocks_per_PSU = 128      # 每次可同時操作的 block 數
start_time = 0
read_request_interval = 31000
write_request_interval = 210000
write_counter = 0
read_counter = 0
MQSim_traces = []
SimpleSSD_traces = []

# ===================== trace 儲存設定 ===================== #
current_dir = Path(__file__).resolve().parent 
ssdconfig_path = current_dir / "ssdconfig.xml"
MQSim_output_dir = current_dir / "traces" / "MQSim.trace"
SimpleSSD_output_dir = current_dir / "traces" / "SimpleSSD.trace"

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
sector_size = 4096
n_sectors_per_page = int(page_capacity / sector_size)

# ===================== MQSim_Trace 類別 ===================== #
class MQSim_Trace:
    def __init__(self, time, level, lsa, size, type):
        self.time = time
        self.level = level
        self.lsa = lsa
        self.size = size
        self.type = type
    def get_str(self):
        return "%d %d %d %d %d" % (self.time, self.level, self.lsa, self.size, self.type)

class SimpleSSD_Trace:
    def __init__(self, time, type, lsa, size):
        self.time = time
        self.lsa = lsa
        self.size = size
        self.type = type
    def get_str(self):
        return "%d %s %d %d" % (self.time, self.type, self.lsa, self.size)

# ===================== 計算寫入與讀取次數 ===================== #
bitvector_length = 4 ** k
bits_per_block = word_lines / 2
blocks_per_bitvector = math.ceil(bitvector_length / bits_per_block)
read_groups = math.ceil(num_reads / reads_per_block)
total_block_required = read_groups * blocks_per_bitvector

if total_block_required >= blocks_per_PSU:
    print("total_block_required >= blocks_per_PSU")
    total_PSU_required = math.ceil(total_block_required / blocks_per_PSU)
    write_block_count = total_block_required
    read_block_count = num_transcripts * total_block_required

else:
    print("total_block_required < blocks_per_PSU")
    total_PSU_required = 1
    write_block_count = total_block_required
    temp = math.ceil(blocks_per_PSU / total_block_required)
    read_block_count = math.ceil(num_transcripts / temp * total_block_required)

# ===================== WRITE =====================
for batch in range(total_PSU_required):
    for channel in range(n_channels): # 這層for以下就是等於1 PSU
        for chip in range(n_chips_per_channel):
            for die in range(n_dies_per_chip):
                if perellelism_level == "channel":
                    chip = 0
                    die = 0
                    batch = 0
                elif perellelism_level == "chip":
                    die = 0
                    batch = 0
                elif perellelism_level == "die":
                    batch = 0
                
                start_lsa = (
                    channel * n_chips_per_channel * n_dies_per_chip * n_planes_per_die * n_blocks_per_plane * n_pages_per_block * n_sectors_per_page # channel offset
                    + chip * n_dies_per_chip * n_planes_per_die * n_blocks_per_plane * n_pages_per_block * n_sectors_per_page # chip offset
                    + die * n_planes_per_die * n_blocks_per_plane * n_pages_per_block * n_sectors_per_page # die offset
                    + batch * n_pages_per_block * n_sectors_per_page # block offset
                )
                
                if write_counter < write_block_count: # 總共要寫total_block_required個block
                    for i in range(n_pages_per_block): # 一次IMS寫入會寫滿一個block = 256個wordline(page)
                        start_lsa += page_capacity * i
                        MQSim_trace = MQSim_Trace(start_time, 0, start_lsa, 1, 0) # Request_Size_In_Sectors = 1，NAND flash write時一樣會write整頁
                        SimpleSSD_trace = SimpleSSD_Trace(start_time, "W", start_lsa, 1)
                        MQSim_traces.append(MQSim_trace)
                        SimpleSSD_traces.append(SimpleSSD_trace)
                        start_time += write_request_interval
                    write_counter += 1
                else:
                    break

# ===================== READ =====================
for batch in range(total_PSU_required * num_transcripts):
    for channel in range(n_channels): # 這層for以下就是等於1 PSU
        for chip in range(n_chips_per_channel):
            for die in range(n_dies_per_chip):
                if perellelism_level == "channel":
                    chip = 0
                    die = 0
                    batch = 0
                elif perellelism_level == "chip":
                    die = 0
                    batch = 0
                elif perellelism_level == "die":
                    batch = 0
                
                start_lsa = (
                    channel * n_chips_per_channel * n_dies_per_chip * n_planes_per_die * n_blocks_per_plane * n_pages_per_block * n_sectors_per_page # channel offset
                    + chip * n_dies_per_chip * n_planes_per_die * n_blocks_per_plane * n_pages_per_block * n_sectors_per_page # chip offset
                    + die * n_planes_per_die * n_blocks_per_plane * n_pages_per_block * n_sectors_per_page # die offset
                    + batch * n_pages_per_block * n_sectors_per_page # block offset
                )
                
                if read_counter < read_block_count: # 總共要搜尋read_block_count個block
                    MQSim_trace = MQSim_Trace(start_time, 0, start_lsa, 1, 1) # Request_Size_In_Sectors = 1，NAND flash read時一樣會read整頁
                    SimpleSSD_trace = SimpleSSD_Trace(start_time, "R", start_lsa, 1)
                    MQSim_traces.append(MQSim_trace)
                    SimpleSSD_traces.append(SimpleSSD_trace)
                    read_counter += 1
                    start_time += read_request_interval
                else:
                    break

print(f"num_reads = {num_reads} , num_transcripts = {num_transcripts}")
print(f"k = {k} , bitvector_length = {bitvector_length} , blocks_per_bitvector = {blocks_per_bitvector}")
print(f"total_block_required = {total_block_required}")
print(f"write_requests = {math.ceil(write_counter * word_lines)}")
print(f"read_requests = {math.ceil(read_counter)}")

# ===================== 輸出 MQSim_Trace =====================
MQSim_output_dir.parent.mkdir(parents=True, exist_ok=True)
with open(MQSim_output_dir, "w") as f:
    for trace in MQSim_traces:
        f.write(trace.get_str() + "\n")

# SimpleSSD_output_dir.parent.mkdir(parents=True, exist_ok=True)
# with open(SimpleSSD_output_dir, "w") as f:
#     for trace in SimpleSSD_traces:
#         f.write(trace.get_str() + "\n")