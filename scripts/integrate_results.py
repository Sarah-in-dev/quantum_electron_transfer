#!/usr/bin/env python3
import json
import os
import csv

print("Starting analysis integration...")

# Check if files exist
quantum_file = "../results/quantum_residues.txt"
conservation_file = "../results/conservation_scores.txt"
distance_file = "../results/cofactor_distances.txt"

# Dictionary to store our results
summary = {}

# Process quantum_residues.txt if it exists
if os.path.exists(quantum_file) and os.path.getsize(quantum_file) > 0:
    print(f"Processing {quantum_file}...")
    
    # Initialize counters
    total_residues = 0
    high_count = 0
    medium_count = 0
    low_count = 0
    high_conservation_sum = 0.0
    overall_conservation_sum = 0.0
    
    # Process file manually
    with open(quantum_file, 'r') as f:
        # Skip header
        next(f)
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:  # Ensure we have enough columns
                total_residues += 1
                
                # Extract conservation score
                try:
                    conservation = float(parts[4])
                    overall_conservation_sum += conservation
                except ValueError:
                    conservation = 0.0
                
                # Check quantum relevance
                relevance = parts[5]
                if relevance == "High":
                    high_count += 1
                    high_conservation_sum += conservation
                elif relevance == "Medium":
                    medium_count += 1
                elif relevance == "Low":
                    low_count += 1
    
    # Calculate means
    mean_conservation_overall = overall_conservation_sum / total_residues if total_residues > 0 else 0
    mean_conservation_high = high_conservation_sum / high_count if high_count > 0 else 0
    
    # Store in summary
    summary["total_residues_analyzed"] = total_residues
    summary["high_relevance_count"] = high_count
    summary["medium_relevance_count"] = medium_count
    summary["low_relevance_count"] = low_count
    summary["mean_conservation_high"] = mean_conservation_high
    summary["mean_conservation_overall"] = mean_conservation_overall
    
    print(f"Processed {total_residues} residues")

# Process cofactor_distances.txt if it exists
if os.path.exists(distance_file) and os.path.getsize(distance_file) > 0:
    print(f"Processing {distance_file}...")
    
    # Initialize variables
    distances = []
    pair_count = 0
    
    # Process file manually
    with open(distance_file, 'r') as f:
        # Skip header
        next(f)
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:  # Ensure we have enough columns
                pair_count += 1
                try:
                    distance = float(parts[2])
                    distances.append(distance)
                except ValueError:
                    pass
    
    # Calculate statistics if we have distances
    if distances:
        min_distance = min(distances)
        max_distance = max(distances)
        mean_distance = sum(distances) / len(distances)
        
        # Store in summary
        summary["cofactor_pairs_analyzed"] = pair_count
        summary["min_cofactor_distance"] = min_distance
        summary["max_cofactor_distance"] = max_distance
        summary["mean_cofactor_distance"] = mean_distance
    
    print(f"Processed {pair_count} cofactor pairs")

# Save summary to JSON
if summary:
    print("Saving summary statistics...")
    with open("../results/summary_statistics.json", "w") as f:
        json.dump(summary, f, indent=2)

# Find top conserved residues near cofactors
if os.path.exists(quantum_file) and os.path.getsize(quantum_file) > 0:
    print("Finding top conserved quantum residues...")
    
    # Read all residues with their conservation scores
    residues = []
    with open(quantum_file, 'r') as f:
        # Skip header
        next(f)
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:  # Ensure we have enough columns
                try:
                    residue_data = {
                        "chain": parts[0],
                        "position": int(parts[1]),
                        "residue": parts[2],
                        "distance_to_cofactor": float(parts[3]),
                        "conservation": float(parts[4]),
                        "quantum_relevance": parts[5]
                    }
                    residues.append(residue_data)
                except (ValueError, IndexError):
                    continue
    
    # Sort by conservation score (highest first)
    residues.sort(key=lambda x: x["conservation"], reverse=True)
    
    # Save top 10 (or fewer if we don't have 10)
    top_count = min(10, len(residues))
    if top_count > 0:
        with open("../results/top_conserved_quantum_residues.csv", "w") as f:
            # Write header
            f.write("chain,position,residue,distance_to_cofactor,conservation,quantum_relevance\n")
            
            # Write data
            for i in range(top_count):
                r = residues[i]
                f.write(f"{r['chain']},{r['position']},{r['residue']},{r['distance_to_cofactor']:.2f},{r['conservation']:.4f},{r['quantum_relevance']}\n")
        
        print(f"Saved top {top_count} conserved residues")

# Print summary to console
print("\nAnalysis complete!")
if "total_residues_analyzed" in summary:
    print(f"Total residues analyzed: {summary['total_residues_analyzed']}")
    high_pct = summary['high_relevance_count']/summary['total_residues_analyzed']*100 if summary['total_residues_analyzed'] > 0 else 0
    med_pct = summary['medium_relevance_count']/summary['total_residues_analyzed']*100 if summary['total_residues_analyzed'] > 0 else 0
    low_pct = summary['low_relevance_count']/summary['total_residues_analyzed']*100 if summary['total_residues_analyzed'] > 0 else 0
    
    print(f"High quantum relevance: {summary['high_relevance_count']} ({high_pct:.1f}%)")
    print(f"Medium quantum relevance: {summary['medium_relevance_count']} ({med_pct:.1f}%)")
    print(f"Low quantum relevance: {summary['low_relevance_count']} ({low_pct:.1f}%)")
    print(f"Mean conservation in high relevance residues: {summary['mean_conservation_high']:.3f}")
    print(f"Overall mean conservation: {summary['mean_conservation_overall']:.3f}")

if "min_cofactor_distance" in summary:
    print(f"Cofactor distances range: {summary['min_cofactor_distance']:.1f} - {summary['max_cofactor_distance']:.1f} Å")
