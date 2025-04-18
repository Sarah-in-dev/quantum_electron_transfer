#!/usr/bin/env python3
from Bio.PDB import PDBParser
import os

print("Analyzing alternative PDB file for cofactors...")

# Check if PDB file exists
pdb_file = "../data/4UB6.pdb"
if not os.path.exists(pdb_file):
    print(f"ERROR: Could not find PDB file {pdb_file}")
    exit(1)

print(f"Using PDB file: {pdb_file}")

# Parse the structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("psii", pdb_file)

# Scan for key cofactors in D1/D2 chains
print("\nScanning for key cofactors in D1 (Chain A) and D2 (Chain D)...")

d1_cofactors = []
d2_cofactors = []

for model in structure:
    for chain in model:
        if chain.id == 'A':  # D1 protein
            for residue in chain:
                if residue.id[0] != ' ':  # Non-standard residue (cofactor)
                    d1_cofactors.append((residue.resname, residue.id[1]))
        
        if chain.id == 'D':  # D2 protein
            for residue in chain:
                if residue.id[0] != ' ':  # Non-standard residue (cofactor)
                    d2_cofactors.append((residue.resname, residue.id[1]))

# Print cofactors
print("\nD1 (Chain A) cofactors:")
for resname, resid in d1_cofactors:
    print(f"  {resname} {resid}")

print("\nD2 (Chain D) cofactors:")
for resname, resid in d2_cofactors:
    print(f"  {resname} {resid}")

# Look specifically for pheophytin in D2
print("\nSearching specifically for pheophytin in D2 (Chain D)...")
for resname, resid in d2_cofactors:
    if resname in ['PHO', 'PHE', 'BPH']:
        print(f"  Found potential pheophytin in D2: {resname} {resid}")
