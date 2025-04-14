#!/usr/bin/env python3
"""
This script generates basic figures for the quantum biology project
using only standard Python libraries and matplotlib.
"""

import os
import matplotlib
# Set non-interactive backend for HPC environments
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

print("Starting figure generation...")

# Check if matplotlib is available
try:
    import matplotlib.pyplot as plt
    has_matplotlib = True
    print("Matplotlib is available - will generate figures")
except ImportError:
    has_matplotlib = False
    print("WARNING: Matplotlib is not available - will skip figure generation")

# Only proceed if matplotlib is available
if has_matplotlib:
    # Create output directory if it doesn't exist
    os.makedirs("../results", exist_ok=True)
    
    # Load quantum residue data if available
    quantum_file = "../results/quantum_residues.txt"
    if os.path.exists(quantum_file) and os.path.getsize(quantum_file) > 0:
        print(f"Loading data from {quantum_file}...")
        
        # Process file manually
        distances = []
        conservations = []
        relevance_categories = []
        
        with open(quantum_file, 'r') as f:
            # Skip header
            next(f)
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:  # Ensure we have enough columns
                    try:
                        distance = float(parts[3])
                        conservation = float(parts[4])
                        relevance = parts[5]
                        
                        distances.append(distance)
                        conservations.append(conservation)
                        relevance_categories.append(relevance)
                    except (ValueError, IndexError):
                        continue
        
        # Count relevance categories
        high_count = relevance_categories.count("High")
        medium_count = relevance_categories.count("Medium")
        low_count = relevance_categories.count("Low")
        
        # Figure 1: Conservation Distribution Histogram
        print("Generating conservation distribution histogram...")
        plt.figure(figsize=(10, 6))
        plt.hist(conservations, bins=20, color='blue', alpha=0.7)
        plt.title("Conservation Score Distribution", fontsize=14)
        plt.xlabel("Conservation Score", fontsize=12)
        plt.ylabel("Count", fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.savefig("../results/figure1_conservation_distribution.png", dpi=300)
        plt.close()
        
        # Figure 2: Scatter plot of distance vs conservation
        print("Generating distance vs conservation scatter plot...")
        plt.figure(figsize=(10, 6))
        
        # Split data by relevance category for color coding
        high_distances = [d for d, r in zip(distances, relevance_categories) if r == "High"]
        high_conservations = [c for c, r in zip(conservations, relevance_categories) if r == "High"]
        
        medium_distances = [d for d, r in zip(distances, relevance_categories) if r == "Medium"]
        medium_conservations = [c for c, r in zip(conservations, relevance_categories) if r == "Medium"]
        
        low_distances = [d for d, r in zip(distances, relevance_categories) if r == "Low"]
        low_conservations = [c for c, r in zip(conservations, relevance_categories) if r == "Low"]
        
        # Plot each category
        plt.scatter(high_distances, high_conservations, color='blue', label='High', alpha=0.7)
        plt.scatter(medium_distances, medium_conservations, color='green', label='Medium', alpha=0.7)
        plt.scatter(low_distances, low_conservations, color='gray', label='Low', alpha=0.7)
        
        plt.title("Relationship Between Cofactor Distance and Conservation", fontsize=14)
        plt.xlabel("Distance to Nearest Cofactor (Å)", fontsize=12)
        plt.ylabel("Conservation Score", fontsize=12)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig("../results/figure2_distance_conservation.png", dpi=300)
        plt.close()
        
        # Figure 3: Bar chart of quantum relevance distribution
        print("Generating quantum relevance distribution bar chart...")
        plt.figure(figsize=(8, 6))
        categories = ['High', 'Medium', 'Low']
        counts = [high_count, medium_count, low_count]
        colors = ['blue', 'green', 'gray']
        
        plt.bar(categories, counts, color=colors, alpha=0.7)
        plt.title("Distribution of Residues by Quantum Relevance Category", fontsize=14)
        plt.xlabel("Quantum Relevance Category", fontsize=12)
        plt.ylabel("Count", fontsize=12)
        
        # Add count labels on top of bars
        for i, count in enumerate(counts):
            plt.text(i, count + 0.5, str(count), ha='center')
            
        plt.tight_layout()
        plt.savefig("../results/figure3_quantum_relevance_distribution.png", dpi=300)
        plt.close()
        
        print("All figures generated successfully!")
    else:
        print(f"ERROR: Could not find or read {quantum_file}")
        print("No figures generated.")

else:
    print("To install matplotlib, you can use:")
    print("  pip install --user matplotlib")
    print("Or setup a virtual environment:")
    print("  python -m venv ~/quantum_env")
    print("  source ~/quantum_env/bin/activate")
    print("  pip install matplotlib")
