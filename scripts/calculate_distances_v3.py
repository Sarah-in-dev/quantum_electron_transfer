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

# Define specific cofactors for D1/D2 reaction center that we know exist
cofactors = {
    # In D1 (Chain A)
    'P680_A': {'chain': 'A', 'resname': 'CLA', 'resid': 405, 'atom': 'MG'},  # Special pair chlorophyll in D1
    'PheoA': {'chain': 'A', 'resname': 'PHO', 'resid': 408, 'atom': 'NA'},   # Pheophytin in D1
    'QA': {'chain': 'A', 'resname': 'PL9', 'resid': 414, 'atom': 'O1'},      # Plastoquinone QA
    
    # In D2 (Chain D)
    'P680_D': {'chain': 'D', 'resname': 'CLA', 'resid': 402, 'atom': 'MG'},  # Special pair chlorophyll in D2
    'QB': {'chain': 'D', 'resname': 'PL9', 'resid': 405, 'atom': 'O1'}       # Plastoquinone QB
}

# Extract additional CLA chlorophylls in D2 that might be part of the electron transfer chain
print("\nSearching for accessory chlorophylls in D2...")
for model in structure:
    for chain in model:
        if chain.id == 'D':
            for residue in chain:
                if residue.resname == 'CLA' and residue.id[1] != 402:  # Not P680_D
                    print(f"Found additional chlorophyll in D2: CLA {residue.id[1]}")
                    cofactors[f'CLA_D_{residue.id[1]}'] = {
                        'chain': 'D', 
                        'resname': 'CLA', 
                        'resid': residue.id[1], 
                        'atom': 'MG'
                    }

# Extract cofactor coordinates
cofactor_coords = {}
for name, info in cofactors.items():
    found = False
    
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

print(f"\nSuccessfully identified {len(cofactor_coords)} cofactors")

# Define primary electron transfer pathways
electron_path = [
    ('P680_A', 'P680_D'),   # Between special pair chlorophylls
    ('P680_A', 'PheoA'),    # Primary electron transfer in D1
    ('PheoA', 'QA'),        # Electron transfer to primary quinone
    ('QA', 'QB')            # Electron transfer between quinones
]

# Add paths for alternative chlorophylls in D2 if found
for name in cofactor_coords:
    if name.startswith('CLA_D_') and name != 'P680_D':
        electron_path.append(('P680_D', name))  # P680_D to accessory CLA
        electron_path.append((name, 'QB'))      # Accessory CLA to QB (alternative path)

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
