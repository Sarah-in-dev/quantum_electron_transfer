## Analysis Pipeline

### 1. Data Collection
- Download D1/D2 protein sequences from UniProt
- Download PSII structures from PDB

### 2. Sequence Analysis
- Multiple sequence alignment with MUSCLE
- Conservation analysis of aligned sequences

### 3. Structural Analysis
- Map conservation onto protein structures
- Calculate distances between electron transfer cofactors
- Identify residues near cofactors

### 4. Quantum Relevance Classification
- Classify residues based on:
  - Distance to cofactors
  - Conservation score
  - Structural position

### 5. Statistical Analysis
- Correlation between conservation and distance
- Comparison of conservation between relevance categories
- Distribution of residues by quantum relevance

## Running the Pipeline
1. Set up project structure: `mkdir -p quantum_project/{data,scripts,results}`
2. Download sequences and structures to `data/` directory
3. Run scripts in the following order:
   - `analyze_conservation.py`
   - `map_conservation.py`
   - `calculate_distances.py`
   - `identify_quantum_residues.py`
   - `statistical_analysis.R`
   - `integrate_results.py`
   - `generate_figures.py`

## Results
The key results files are:
- `quantum_residues.txt`: Classification of residues by quantum relevance
- `summary_statistics.json`: Overview of key findings
- Figure files for publication

## References
- Engel GS, et al. (2007). Evidence for wavelike energy transfer through quantum coherence in photosynthetic systems. Nature, 446(7137), 782-786.
- Romero E, et al. (2014). Quantum coherence in photosynthesis for efficient solar-energy conversion. Nature Physics, 10(9), 676-682.
- Umena Y, et al. (2011). Crystal structure of oxygen-evolving photosystem II at a resolution of 1.9 Å. Nature, 473(7345), 55-60.
EOF
