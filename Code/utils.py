import os
import chardet
import xml.etree.ElementTree as ET
import xml.dom.minidom
import numpy as np
from bidict import bidict
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict, Counter

num_read_dict = bidict({
    "0.1M": 100_000,
    "1M": 1_000_000,
    "3.2M": 3_200_000,
    "10M": 10_000_000,
    "100M": 100_000_000
})

# ===================== Trace 類別 ===================== #
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

def generate_workload_xml(num_traces, output_path, ssdconfig_path=None):
    # 如果沒傳就用預設路徑（utils.py 上一層的 SSD_config.xml）
    if ssdconfig_path is None:
        ssdconfig_path = Path(__file__).resolve().parent.parent / "SSD_config.xml"

    # 讀取 SSD_config.xml
    tree = ET.parse(ssdconfig_path)
    root = tree.getroot()
    device_para_node = root.find("Device_Parameter_Set")
    flash_para_node = device_para_node.find("Flash_Parameter_Set")

    n_channels = int(device_para_node.find("Flash_Channel_Count").text)
    n_chips_per_channel = int(device_para_node.find("Chip_No_Per_Channel").text)
    n_dies_per_chip = int(flash_para_node.find("Die_No_Per_Chip").text)
    n_planes_per_die = int(flash_para_node.find("Plane_No_Per_Die").text)

    # 生成 ID 範圍字串
    channel_ids = ",".join(str(i) for i in range(n_channels))
    chip_ids = ",".join(str(i) for i in range(n_chips_per_channel))
    die_ids = ",".join(str(i) for i in range(n_dies_per_chip))
    plane_ids = ",".join(str(i) for i in range(n_planes_per_die))

    # 建立 XML 結構
    root_elem = ET.Element("MQSim_IO_Scenarios")

    for i in range(num_traces):
        scenario = ET.SubElement(root_elem, "IO_Scenario")
        flow = ET.SubElement(scenario, "IO_Flow_Parameter_Set_Trace_Based")

        ET.SubElement(flow, "Priority_Class").text = "HIGH"
        ET.SubElement(flow, "Device_Level_Data_Caching_Mode").text = "WRITE_READ_CACHE"
        ET.SubElement(flow, "Channel_IDs").text = channel_ids
        ET.SubElement(flow, "Chip_IDs").text = chip_ids
        ET.SubElement(flow, "Die_IDs").text = die_ids
        ET.SubElement(flow, "Plane_IDs").text = plane_ids
        ET.SubElement(flow, "Initial_Occupancy_Percentage").text = "0"
        ET.SubElement(flow, "File_Path").text = f"traces/MQSim_{i}.trace"
        ET.SubElement(flow, "Percentage_To_Be_Executed").text = "100"
        ET.SubElement(flow, "Relay_Count").text = "1"
        ET.SubElement(flow, "Time_Unit").text = "NANOSECOND"

    # 輸出漂亮格式的 XML
    rough_string = ET.tostring(root_elem, encoding="us-ascii")
    reparsed = xml.dom.minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="\t")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="us-ascii") as f:
        f.write(pretty_xml)

# def generate_workload_xml(num_traces, output_path):
#     root = ET.Element("MQSim_IO_Scenarios")

#     for i in range(num_traces):
#         scenario = ET.SubElement(root, "IO_Scenario")
#         flow = ET.SubElement(scenario, "IO_Flow_Parameter_Set_Trace_Based")

#         ET.SubElement(flow, "Priority_Class").text = "HIGH"
#         ET.SubElement(flow, "Device_Level_Data_Caching_Mode").text = "WRITE_READ_CACHE"
#         ET.SubElement(flow, "Channel_IDs").text = "0,1,2,3,4,5,6,7"
#         ET.SubElement(flow, "Chip_IDs").text = "0,1,2,3,4,5,6,7"
#         ET.SubElement(flow, "Die_IDs").text = "0,1"
#         ET.SubElement(flow, "Plane_IDs").text = "0"
#         ET.SubElement(flow, "Initial_Occupancy_Percentage").text = "0"
#         ET.SubElement(flow, "File_Path").text = f"traces/MQSim_{i}.trace"
#         ET.SubElement(flow, "Percentage_To_Be_Executed").text = "100"
#         ET.SubElement(flow, "Relay_Count").text = "1"
#         ET.SubElement(flow, "Time_Unit").text = "NANOSECOND"

#     rough_string = ET.tostring(root, encoding="us-ascii")
#     reparsed = xml.dom.minidom.parseString(rough_string)
#     pretty_xml = reparsed.toprettyxml(indent="\t")

#     output_path.parent.mkdir(parents=True, exist_ok=True)
#     with open(output_path, "w", encoding="us-ascii") as f:
#         f.write(pretty_xml)

def count_reads(folder_path):
    folder = Path(folder_path)
    fastq_files = list(folder.glob("*.fastq")) + list(folder.glob("*.fq"))
    total_reads = 0
    for file_path in fastq_files:
        with open(file_path, 'r') as f:
            total_reads += sum(1 for _ in f) // 4
    return total_reads

def count_transcripts(folder_path):
    folder = Path(folder_path)
    fasta_files = list(folder.glob("*.fasta")) + list(folder.glob("*.fa"))
    total = 0
    for f in fasta_files:
        # 嘗試自動偵測檔案編碼
        with open(f, 'rb') as raw_file:
            raw_data = raw_file.read(10000)
            result = chardet.detect(raw_data)
            encoding = result['encoding']
        
        # 用偵測到的編碼開啟檔案
        with open(f, 'r', encoding=encoding) as file:
            total += sum(1 for line in file if line.startswith('>'))
    return total

def output_trace(output_dir, traces):
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with open(output_dir, "w") as f:
        for trace in traces:
            f.write(trace.get_str() + "\n")

def load_transcripts(folder_path):
    print("[STEP] Load transcripts")
    transcripts = []
    for file in os.listdir(folder_path):
        if file.endswith(".fasta") or file.endswith(".fa"):
            file_path = os.path.join(folder_path, file)
            for record in SeqIO.parse(file_path, "fasta"):
                transcripts.append(str(record.seq).upper())
    return transcripts

def load_reads(folder_path, k, max_reads=None):
    print("[STEP] Load reads")
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

def extract_all_kmers(transcripts, reads, k):
    kmer_counter = Counter()
    for seq in transcripts + reads:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if set(kmer).issubset({'A', 'C', 'G', 'T'}):
                kmer_counter[kmer] += 1
    return sorted(kmer_counter.keys())