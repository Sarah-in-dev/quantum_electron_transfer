#!/usr/bin/env python3
from Bio.PDB import PDBParser, Selection, NeighborSearch
import numpy as np

# Load structure
parser = PDBParser()
structure = parser.get_structure("psii", "../data/3WU2.pdb")

# Extract all atoms
atoms = list(structure.get_atoms())

# Create neighbor search object
ns = NeighborSearch(atoms)

# Define cofactors
cofactors = []
cofactor_names = {
    'CLA': 'Chlorophyll',
    'PHO': 'Pheophytin',
    'PL9': 'Plastoquinone'
}

for model in structure:
    for chain in model:
        if chain.id in ['A', 'D']:  # D1 and D2 chains
            for residue in chain:
                if residue.resname in cofactor_names:
                    for atom in residue:
                        cofactors.append((atom, residue.resname, residue.id[1], chain.id))

print(f"Found {len(cofactors)} cofactor atoms")

# Find residues near cofactors and calculate distances
quantum_residues = []
for atom, resname, resid, chain in cofactors:
    # Find all atoms within 5Å of this cofactor atom
    neighbors = ns.search(atom.coord, 5.0)
    for neighbor in neighbors:
        res = neighbor.get_parent()
        if res.resname not in cofactor_names and res.id[0] == ' ':  # Standard amino acid
            # Calculate distance between this residue and the cofactor
            distance = np.linalg.norm(neighbor.coord - atom.coord)
            quantum_residues.append((res.get_parent().id, res.id[1], res.resname, distance))

print(f"Found {len(quantum_residues)} residues near cofactors")

# Remove duplicates (keep minimum distance)
unique_residues = {}
for chain, pos, resname, distance in quantum_residues:
    key = (chain, pos, resname)
    if key not in unique_residues or distance < unique_residues[key]:
        unique_residues[key] = distance

print(f"After removing duplicates, found {len(unique_residues)} unique residues")

# Load conservation scores
conservation = {}
try:
    with open("../results/conservation_scores.txt") as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                pos, score = parts
                conservation[int(pos)] = float(score)
    print(f"Loaded {len(conservation)} conservation scores")
except Exception as e:
    print(f"Error loading conservation scores: {e}")
    # Create dummy conservation scores if file not available
    for chain, pos, _ in unique_residues.keys():
        conservation[pos] = 0.5

# Classify residues based on quantum relevance
with open("../results/quantum_residues.txt", "w") as f:
    f.write("chain\tposition\tresidue\tdistance_to_cofactor\tconservation\tquantum_relevance\n")
    for (chain, pos, resname), distance in unique_residues.items():
        # Get conservation score with default of 0.5 if not found
        cons_score = conservation.get(pos, 0.5)
        
        # Classify quantum relevance
        if distance < 14 and cons_score > 0.8:
            relevance = "High"
        elif distance < 20 and cons_score > 0.6:
            relevance = "Medium"
        else:
            relevance = "Low"
        
        f.write(f"{chain}\t{pos}\t{resname}\t{distance:.2f}\t{cons_score:.4f}\t{relevance}\n")

print("Analysis complete - results written to ../results/quantum_residues.txt")
