#!/usr/bin/env python3
from Bio import AlignIO
from Bio.Align import AlignInfo
import numpy as np

# Load the alignment
alignment = AlignIO.read("../results/aligned_proteins.fasta", "fasta")

# Calculate conservation at each position
summary_align = AlignInfo.SummaryInfo(alignment)
conservation_dict = {}
for i in range(alignment.get_alignment_length()):
    column = alignment[:, i]
    # Count frequency of most common amino acid
    aa_count = {}
    for aa in column:
        if aa not in aa_count:
            aa_count[aa] = 0
        aa_count[aa] += 1
    
    max_freq = max(aa_count.values()) / len(column)
    conservation_dict[i] = max_freq

# Write conservation scores to file
with open("../results/conservation_scores.txt", "w") as f:
    f.write("position\tconservation_score\n")
    for pos, score in conservation_dict.items():
        f.write(f"{pos+1}\t{score:.4f}\n")

# Identify highly conserved regions (score > 0.8)
conserved_regions = []
current_region = []
for pos, score in conservation_dict.items():
    if score > 0.8:
        current_region.append(pos)
    elif current_region:
        if len(current_region) > 3:  # Minimum region length
            conserved_regions.append((min(current_region), max(current_region)))
        current_region = []

# Write conserved regions to file
with open("../results/conserved_regions.txt", "w") as f:
    f.write("start\tend\tlength\n")
    for start, end in conserved_regions:
        f.write(f"{start+1}\t{end+1}\t{end-start+1}\n")

