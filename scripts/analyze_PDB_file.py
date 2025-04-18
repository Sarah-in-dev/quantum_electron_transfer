#!/usr/bin/env python3
from Bio.PDB import PDBParser
import os

print("Analyzing PDB structure to identify cofactors...")

# Check if PDB file exists
pdb_file = "../data/3WU2.pdb"
if not os.path.exists(pdb_file):
    pdb_file = "../data/4UB6.pdb"
    if not os.path.exists(pdb_file):
        print("ERROR: Could not find PDB files")
        exit(1)

print(f"Using PDB file: {pdb_file}")

# Parse the structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("psii", pdb_file)

# Scan for hetero residues (potential cofactors)
het_residues = {}

for model in structure:
    print(f"Model: {model.id}")
    for chain in model:
        chain_het = []
        for residue in chain:
            # Hetero residues are non-amino acid molecules
            if residue.id[0] != ' ':
                res_info = {
                    'resname': residue.resname,
                    'resid': residue.id[1],
                    'atoms': len(list(residue.get_atoms())),
                    'atom_names': [atom.name for atom in residue][:10]  # First 10 atoms
                }
                chain_het.append(res_info)
        
        if chain_het:
            het_residues[chain.id] = chain_het
            print(f"  Chain {chain.id}: Found {len(chain_het)} hetero residues")

# Print details of potential cofactors
print("\nPotential cofactors:")
for chain_id, residues in het_residues.items():
    print(f"\nChain {chain_id}:")
    for res in residues:
        # Focus on larger molecules (likely cofactors)
        if res['atoms'] > 10:
            print(f"  {res['resname']} {res['resid']} - {res['atoms']} atoms")
            print(f"    Sample atoms: {', '.join(res['atom_names'])}")
            
            # Identify likely chlorophylls (have MG), pheophytins, and quinones
            if 'MG' in res['atom_names']:
                print(f"    *** Likely CHLOROPHYLL (has MG atom)")
            elif res['atoms'] > 30 and res['resname'] in ['PHO', 'BPH', 'PHE']:
                print(f"    *** Likely PHEOPHYTIN")
            elif res['atoms'] < 30 and any(name in res['atom_names'] for name in ['O1', 'O2']):
                print(f"    *** Likely QUINONE")
