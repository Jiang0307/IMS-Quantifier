import pandas as pd
import numpy as np
from numpy import genfromtxt
transcriptome = pd.read_csv('C:\Users\User\Desktop\GitHub\IMS-quantifier\Dataset\Transcript\transcirpts.fasta', header=None, delimiter=',')
RNA_read = pd.read_csv("C:\Users\User\Desktop\GitHub\IMS-quantifier\Dataset\Read\SRR390728_1.fastq", header=None)

# Segment name for reads
read_seq_name = np.array(RNA_read.iloc[::2,:])
for i in range(len(read_seq_name)):
    first = read_seq_name[i][0].split("/")
    second = first[1].split(";")[0]
    read_seq_name[i][0] = second[0:16]

# Read Sequence
read_seq = np.array(read.iloc[1::2,:])

# Turn it into dataframe
RNA_read_df = pd.DataFrame(read_seq_name, columns=["Reference Region Name"])
RNA_read_df["Read Sequence"] = read_seq

# Transcriptome dict to refine fasta file
transcriptome_dict = {}
storage = ""

for i in range(len(transcriptome.values)):
    # Reference Segment Name
    if transcriptome.values[i][0][0] == ">":
        storage = transcriptome.values[i][0][1:]
        if storage not in transcriptome_dict.keys():
            transcriptome_dict[storage] = ""
    # Sequence
    else: 
        transcriptome_dict[storage] += transcriptome.values[i][0]

# Turn it into Dataframe
transcriptome_df = pd.DataFrame(zip(list(transcriptome_dict.keys()), list(transcriptome_dict.values())),
                        columns=["Reference Name", "Sequence"])

# transcriptome_df_rev = transcriptome_df.copy()
# for i in range(len(transcriptome_df_rev)):
#     transcriptome_df_rev.iloc[i][1] = transcriptome_df_rev.iloc[i][1][::-1]
    
RNA_read_df_rev = RNA_read_df.copy()    
for i in range(len(RNA_read_df_rev)):
    RNA_read_df_rev.iloc[i][1] = RNA_read_df_rev.iloc[i][1][::-1]
    target_sequence = list(RNA_read_df_rev.iloc[i][1])
    for j in range(len(target_sequence)):
        if target_sequence[j] == "C":
            target_sequence[j] = "G"
            
        elif target_sequence[j] == "A":
            target_sequence[j] = "T"
            
        elif target_sequence[j] == "G":
            target_sequence[j] = "C"
            
        elif target_sequence[j] == "T":
            target_sequence[j] = "A"
    RNA_read_df_rev.iloc[i][1] = "".join(target_sequence)


def kmer_composition(seq, k):
    """
    This function finds Kmer composition length of k from the given sequence
    """
    all_kmers = list()
    for i in range(len(seq) - k + 1):
        current_kmer = seq[i:(i + k)]
        all_kmers.append(current_kmer)
    return all_kmers

def common_elements(combined_list):
    """
    This function finds common elements in the list and return list with unique elements
    """
    if len(combined_list) == 0:
        return []
    
    unique = set(combined_list[0])
    for l in combined_list[1:]:
        unique = unique.intersection(l)

    return list(unique)

def Hash(reads, ref, k):
    """
    This function returns the hash table for transcriptome. 
    Each key is the kmer, and its corresponding value is the reference name from which
    the kmer came from.
    """
    # Store edges and reference
    edges = {}
    reference = {}
    
    # To keep track of reads being examined
    nth_read = 0
    
    # Loop through each sequence 
    for read in reads:
        i = 0
        while i+k < len(read):
            # left kmer
            left_mer = read[i:i+k]
            # right kmer
            right_mer = read[i+1:i+k+1]
            
            # If already in the edges
                # This indicates new sequence or circular graph formation
            if left_mer in edges.keys():
                # If reference is in the sequence: same node from difference segment
                if ref[nth_read] not in reference[left_mer]:
                    reference[left_mer] += [ref[nth_read]]
                edges[left_mer] += [left_mer + right_mer[-1]]
                    
            # If node does not exist, start of a new sequecne
            else:
                edges[left_mer] = [left_mer + right_mer[-1]]
                reference[left_mer] = [ref[nth_read]]
            
            if right_mer not in edges.keys():
                edges[right_mer] = []
                reference[right_mer] = [ref[nth_read]]

            i += 1
        nth_read += 1

    return reference
        
def pseudoalignment(reference, reads, k):
    """
    This function aligns the read to the given hash table (reference).
    """
    read_aligned = {}
    nth_read = 0
    
    # loops through each read
    for read in reads:
        # used as key to keep track of reads
        key = "read" + str(nth_read)
        read_aligned[key] = []
        
        # Find kmer composition
        compositions = kmer_composition(read, k)
        
        # For each kmer, try to align 
        for kmer in compositions:
            # If cannot be aligned, break out of the loop
            if kmer not in reference.keys():
                # clear it
                read_aligned[key] = []
                break
            # If there is a match, add its reference to the value in the dictionary
            else:
                read_aligned[key] += [reference[kmer]]
                
        # Find common elements of equivalence class
        read_aligned[key] = common_elements(read_aligned[key])
        nth_read += 1

    return read_aligned



final_table



# Get Hash Table from the transcriptome
hash_table = Hash(transcriptome_df["Sequence"][:], transcriptome_df["Reference Name"][:], 31)

# Get alignment results for forward and reverse strand of reads
aligned_result_for = pseudoalignment(hash_table, RNA_read_df["Read Sequence"][:], 31)
aligned_result_rev = pseudoalignment(hash_table, RNA_read_df_rev["Read Sequence"][:], 31)

# Combine alignment results for forward and reverse alignment
combined = {}
for key in aligned_result_for.keys():
    forward_reference = aligned_result_for[key]
    reverse_reference = aligned_result_rev[key]
    final_result = set(forward_reference + reverse_reference)
    combined[key] = list(final_result)
    
# Listed equivalence to string
for key in combined.keys():
    combined[key] = ', '.join(combined[key])
    
equiv_list = list(combined.values())
set_equiv_list = list(set(equiv_list))

# Count number of reads per equivalence class
equiv_read_count = {}
for i in range(len(set_equiv_list)):
    equiv_read_count[set_equiv_list[i]] = equiv_list.count(set_equiv_list[i])
    
equivalence_class = list(equiv_read_count.keys())
read_count = list(equiv_read_count.values())

# number of isomers in each equivalence class
n_isomers = []
for i in range(len(equivalence_class)):
    if len(equivalence_class[i]) == 0:
        n_isomers.append(0)
    else:
        n_isomers.append(equivalence_class[i].count(',') + 1)

# Make a table
final_table = pd.DataFrame(zip(read_count, n_isomers, equivalence_class), columns=["Read Count", "# Isoforms", "Equivalence Class"])



final_table

from matplotlib import pyplot as plt
import seaborn as sns

not_aligned = [(final_table["Read Count"][0])]
aligned = [len(RNA_read_df) - final_table["Read Count"][0]]
column = ["Not Aligned", "Aligned"]
comparison = pd.DataFrame(zip(not_aligned + aligned, column), columns=["Alignment", "Type of Reads"])

# Graph
fgrid = sns.displot(data=final_table, x="# Isoforms")
fgrid2 = sns.catplot(data=comparison, x="Type of Reads", y="Alignment", kind="bar")

# Title
fgrid.fig.suptitle('Distribution of Isoforms in Equivalence Class')
fgrid2.fig.suptitle('% Aligned = 81.93 %')


aligned_percentage = (final_table["Read Count"][1:].sum() / final_table["Read Count"][:].sum()) * 100
aligned_percentage