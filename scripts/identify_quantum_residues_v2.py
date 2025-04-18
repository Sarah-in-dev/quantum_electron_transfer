#!/usr/bin/env python3
from Bio.PDB import PDBParser, Selection, NeighborSearch
import numpy as np
import os

# Create results directory if it doesn't exist
os.makedirs("../results", exist_ok=True)

print("Starting quantum residue identification...")

# Check if PDB file exists
pdb_file = "../data/3WU2.pdb"
if not os.path.exists(pdb_file):
    print(f"ERROR: PDB file not found at {pdb_file}")
    exit(1)

print(f"Using PDB file: {pdb_file}")

# Define known cofactors based on our distance calculation
cofactors = {
    # In D1 (Chain A)
    'P680_A': {'chain': 'A', 'resname': 'CLA', 'resid': 405},
    'PheoA': {'chain': 'A', 'resname': 'PHO', 'resid': 408},
    'QA': {'chain': 'A', 'resname': 'PL9', 'resid': 414},
    
    # In D2 (Chain D)
    'P680_D': {'chain': 'D', 'resname': 'CLA', 'resid': 402},
    'CLA_D_403': {'chain': 'D', 'resname': 'CLA', 'resid': 403},
    'QB': {'chain': 'D', 'resname': 'PL9', 'resid': 405}
}

# Parse the structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("psii", pdb_file)

# Extract all atoms
atoms = list(structure.get_atoms())
print(f"Extracted {len(atoms)} atoms from structure")

# Create neighbor search object
ns = NeighborSearch(atoms)

# Extract cofactor residues
cofactor_residues = []
for name, info in cofactors.items():
    found = False
    for model in structure:
        for chain in model:
            if chain.id == info['chain']:
                for residue in chain:
                    if (residue.resname == info['resname'] and 
                        residue.id[1] == info['resid']):
                        cofactor_residues.append((residue, name))
                        found = True
                        print(f"Found cofactor {name}: {info['resname']} {info['resid']} in chain {info['chain']}")
    
    if not found:
        print(f"WARNING: Could not find cofactor {name}")

print(f"Found {len(cofactor_residues)} cofactor residues")

# Find residues near cofactors and calculate distances
quantum_residues = []
for cofactor, name in cofactor_residues:
    # Use all atoms in the cofactor to search for nearby residues
    for cofactor_atom in cofactor:
        # Find all atoms within 14Å of this cofactor atom (range for quantum effects)
        neighbors = ns.search(cofactor_atom.coord, 14.0)
        for neighbor in neighbors:
            res = neighbor.get_parent()
            # Only consider standard amino acids from D1/D2 chains
            if (res.id[0] == ' ' and 
                res.get_parent().id in ['A', 'D'] and
                res.resname not in ['CLA', 'PHO', 'PL9']):
                
                # Calculate distance
                distance = np.linalg.norm(neighbor.coord - cofactor_atom.coord)
                quantum_residues.append((
                    res.get_parent().id,  # Chain
                    res.id[1],            # Position
                    res.resname,          # Residue type
                    distance,             # Distance
                    name                  # Cofactor name
                ))

print(f"Found {len(quantum_residues)} residue-cofactor interactions")

# Remove duplicates (keep minimum distance)
unique_residues = {}
for chain, pos, resname, distance, cofactor in quantum_residues:
    key = (chain, pos, resname)
    if key not in unique_residues or distance < unique_residues[key][0]:
        unique_residues[key] = (distance, cofactor)

print(f"After removing duplicates, found {len(unique_residues)} unique residues")

# Load conservation scores
conservation = {}
conservation_file = "../results/conservation_scores.txt"
if os.path.exists(conservation_file) and os.path.getsize(conservation_file) > 0:
    try:
        with open(conservation_file) as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    pos, score = parts
                    conservation[int(pos)] = float(score)
        print(f"Loaded {len(conservation)} conservation scores")
    except Exception as e:
        print(f"Error loading conservation scores: {e}")

# Classify residues based on quantum relevance
with open("../results/quantum_residues.txt", "w") as f:
    f.write("chain\tposition\tresidue\tdistance_to_cofactor\tconservation\tquantum_relevance\tcofactor\n")
    for (chain, pos, resname), (distance, cofactor) in unique_residues.items():
        # Get conservation score with default of 0.5 if not found
        cons_score = conservation.get(pos, 0.5)
        
        # Classify quantum relevance based on distance and conservation
        # Adjust thresholds based on the actual cofactor distances we found
        if distance < 8.0 and cons_score > 0.8:
            relevance = "High"
        elif distance < 14.0 and cons_score > 0.6:
            relevance = "Medium"
        else:
            relevance = "Low"
        
        f.write(f"{chain}\t{pos}\t{resname}\t{distance:.2f}\t{cons_score:.4f}\t{relevance}\t{cofactor}\n")

print("Quantum residue identification complete - results written to ../results/quantum_residues.txt")
