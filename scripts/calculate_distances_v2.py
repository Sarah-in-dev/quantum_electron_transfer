#!/usr/bin/env python3
from Bio.PDB import PDBParser
import numpy as np
import os

# Create results directory if it doesn't exist
os.makedirs("../results", exist_ok=True)

print("Calculating distances between electron transfer cofactors...")

# Define the PDB file to use
pdb_file = "../data/3WU2.pdb"
if not os.path.exists(pdb_file):
    print(f"ERROR: Could not find PDB file {pdb_file}")
    exit(1)

print(f"Using PDB file: {pdb_file}")

# Parse the structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("psii", pdb_file)

# Define specific cofactors for D1/D2 reaction center
# Based on PDB analysis
cofactors = {
    # In D1 (Chain A)
    'P680_A': {'chain': 'A', 'resname': 'CLA', 'resid': 405, 'atom': 'MG'},  # Special pair chlorophyll in D1
    'PheoA': {'chain': 'A', 'resname': 'PHO', 'resid': 408, 'atom': 'NA'},   # Pheophytin in D1 (no MG, use NA)
    'QA': {'chain': 'A', 'resname': 'PL9', 'resid': 414, 'atom': 'O1'},      # Plastoquinone QA
    
    # In D2 (Chain D)
    'P680_D': {'chain': 'D', 'resname': 'CLA', 'resid': 402, 'atom': 'MG'},  # Special pair chlorophyll in D2
    'PheoD': {'chain': 'D', 'resname': 'PHO', 'resid': 0, 'atom': 'NA'},     # Pheophytin in D2 (ID to be found)
    'QB': {'chain': 'D', 'resname': 'PL9', 'resid': 405, 'atom': 'O1'}       # Plastoquinone QB
}

# First, let's find the pheophytin in Chain D if it exists
# (it wasn't clearly identified in the PDB analysis)
for model in structure:
    for chain in model:
        if chain.id == 'D':
            for residue in chain:
                if residue.resname == 'PHO':
                    pheo_id = residue.id[1]
                    cofactors['PheoD']['resid'] = pheo_id
                    print(f"Found PheoD: PHO {pheo_id} in Chain D")
                    break

# Extract cofactor coordinates
cofactor_coords = {}
for name, info in cofactors.items():
    found = False
    
    # Skip if the residue ID is 0 (not found)
    if info['resid'] == 0:
        print(f"Skipping {name}: No residue ID defined")
        continue
        
    for model in structure:
        for chain in model:
            if chain.id == info['chain']:
                for residue in chain:
                    if (residue.resname == info['resname'] and 
                        residue.id[1] == info['resid']):
                        
                        # First try to find the specified atom
                        for atom in residue:
                            if atom.name == info['atom']:
                                cofactor_coords[name] = atom.coord
                                found = True
                                print(f"Found {name} at position {atom.coord}")
                                break
                        
                        # If not found (especially for pheophytin which has no MG)
                        # try alternative central atoms
                        if not found and info['resname'] == 'PHO':
                            for atom_name in ['NA', 'NB', 'NC', 'ND', 'C1', 'C2']:
                                for atom in residue:
                                    if atom.name == atom_name:
                                        cofactor_coords[name] = atom.coord
                                        found = True
                                        print(f"Found {name} (using {atom_name}) at position {atom.coord}")
                                        break
                                if found:
                                    break
                        
                        # If still not found, use the first atom
                        if not found:
                            atoms = list(residue.get_atoms())
                            if atoms:
                                cofactor_coords[name] = atoms[0].coord
                                found = True
                                print(f"Found {name} (using first atom) at position {atoms[0].coord}")

print(f"\nSuccessfully identified {len(cofactor_coords)} of {len(cofactors)} cofactors")

# Define electron transfer pathways
electron_path = [
    ('P680_A', 'P680_D'),   # Between special pair chlorophylls
    ('P680_A', 'PheoA'),    # Primary electron transfer in D1
    ('P680_D', 'PheoD'),    # Alternative electron transfer in D2
    ('PheoA', 'QA'),        # Electron transfer to primary quinone
    ('PheoD', 'QB'),        # Electron transfer to secondary quinone
    ('QA', 'QB')            # Electron transfer between quinones
]

# Calculate distances between cofactors in electron transfer paths
with open("../results/cofactor_distances.txt", "w") as f:
    f.write("cofactor1\tcofactor2\tdistance\n")
    
    for pair in electron_path:
        if pair[0] in cofactor_coords and pair[1] in cofactor_coords:
            dist = np.linalg.norm(cofactor_coords[pair[0]] - cofactor_coords[pair[1]])
            f.write(f"{pair[0]}\t{pair[1]}\t{dist:.2f}\n")
            print(f"Distance between {pair[0]} and {pair[1]}: {dist:.2f} Å")
        else:
            missing = []
            if pair[0] not in cofactor_coords:
                missing.append(pair[0])
            if pair[1] not in cofactor_coords:
                missing.append(pair[1])
            print(f"Cannot calculate distance between {pair[0]} and {pair[1]} - missing: {', '.join(missing)}")
            f.write(f"{pair[0]}\t{pair[1]}\tNA\n")

print("\nCofactor distance analysis complete.")
