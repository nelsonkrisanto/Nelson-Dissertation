#!/bin/bash
#SBATCH --job-name=3.6.split_fasta
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10

# Directory where the Python script is located
script_dir="/home/people/23203786/scratch/Nelson-Dissertation/scripts"
script_file="split_fasta.py"

# Directory where the FASTA file is located
fasta_dir="/home/people/23203786/scratch/Nelson-Dissertation/fasta"
fasta_file="1.1.dengue_virus_sequences.fasta"

# Output directory for the split FASTA files
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/fasta"

# Load Python module
module load python/3.9.15

# Change to the script directory
cd "$script_dir"

# Run the Python script
python "$script_file"

# Unload Python module
module unload python/3.9.15