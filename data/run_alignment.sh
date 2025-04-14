#!/bin/bash
#SBATCH --job-name=muscle_align
#SBATCH --output=muscle_align_%j.out
#SBATCH --error=muscle_align_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

module load muscle
muscle -in data/photosystem_proteins.fasta -out results/aligned_proteins.fasta
EOF
