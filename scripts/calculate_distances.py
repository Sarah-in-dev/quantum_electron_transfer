#!/usr/bin/env python3
from Bio.PDB import PDBParser
import numpy as np

# Define cofactor names and their atoms for distance calculation
cofactors = {
    'P680_A': {'chain': 'A', 'resname': 'CLA', 'resid': 1001, 'atom': 'MG'},
    'P680_D': {'chain': 'D', 'resname': 'CLA', 'resid': 1001, 'atom': 'MG'},
    'PheoA': {'chain': 'A', 'resname': 'PHO', 'resid': 1002, 'atom': 'MG'},
    'PheoD': {'chain': 'D', 'resname': 'PHO', 'resid': 1002, 'atom': 'MG'},
    'QA': {'chain': 'A', 'resname': 'PL9', 'resid': 1003, 'atom': 'O1'},
    'QB': {'chain': 'D', 'resname': 'PL9', 'resid': 1004, 'atom': 'O1'},
}

# Load structure
parser = PDBParser()
structure = parser.get_structure("psii", "../data/3WU2.pdb")

# Extract cofactor coordinates
cofactor_coords = {}
for name, info in cofactors.items():
    for model in structure:
        for chain in model:
            if chain.id == info['chain']:
                for residue in chain:
                    if residue.resname == info['resname'] and residue.id[1] == info['resid']:
                        for atom in residue:
                            if atom.name == info['atom']:
                                cofactor_coords[name] = atom.coord

# Calculate distances between cofactors
distances = {}
electron_path = [
    ('P680_A', 'P680_D'),
    ('P680_A', 'PheoA'),
    ('P680_D', 'PheoD'),
    ('PheoA', 'QA'),
    ('PheoD', 'QB'),
    ('QA', 'QB')
]

with open("../results/cofactor_distances.txt", "w") as f:
    f.write("cofactor1\tcofactor2\tdistance\n")
    for pair in electron_path:
        if pair[0] in cofactor_coords and pair[1] in cofactor_coords:
            dist = np.linalg.norm(cofactor_coords[pair[0]] - cofactor_coords[pair[1]])
            distances[pair] = dist
            f.write(f"{pair[0]}\t{pair[1]}\t{dist:.2f}\n")
