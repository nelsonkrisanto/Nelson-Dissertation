#!/bin/bash
#SBATCH --job-name=combine_fasta  # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie,nelson_krisanto@hotmail.com # Where to send mail
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt  # Error log
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt  # Standard output log
#SBATCH --cpus-per-task=5

# Define paths
my_dir="/home/people/23203786/scratch/Nelson-Dissertation"
orthoflavivirus_fasta="$my_dir/raw_data/Orthoflavivirus/orthoflavivirus.fasta"
clustered_fasta="$my_dir/results/fasta/clustered_sequences.fasta"
combined_fasta="$my_dir/results/fasta/combined_sequences.fasta"

# Combine the two FASTA files
cat "$orthoflavivirus_fasta" "$clustered_fasta" > "$combined_fasta"

echo "Combined FASTA file created at: $combined_fasta"
