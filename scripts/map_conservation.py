#!/usr/bin/env python3
import os
from Bio.PDB import PDBParser, Selection, PDBIO
import numpy as np

# Load conservation scores
conservation = {}
with open("../results/conservation_scores.txt") as f:
    next(f)  # Skip header
    for line in f:
        pos, score = line.strip().split()
        conservation[int(pos)] = float(score)

# Load structure
parser = PDBParser()
structure = parser.get_structure("psii", "../data/3WU2.pdb")

# Map conservation scores to B-factor field for visualization
for model in structure:
    for chain in model:
        if chain.id in ['A', 'D']:  # D1 and D2 chains
            for residue in chain:
                if residue.id[0] == ' ':  # Standard amino acid
                    res_num = residue.id[1]
                    if res_num in conservation:
                        for atom in residue:
                            atom.bfactor = conservation[res_num] * 100

# Save modified structure for visualization
io = PDBIO()
io.set_structure(structure)
io.save("../results/conservation_mapped.pdb")
