#!/bin/bash
#SBATCH --job-name=3.6.3.semi_nested_u
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Directory where the TSV files are located
tsv_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv"

# Load Python module
module load python/3.9.15

# Change to the directory where the TSV files are located
cd "$tsv_dir"

# Path to the semi-nested primer combinations script (update this path to your actual semi-nested PCR script location)
semi_nested_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/semi_nested_u.py"

# Define amplicon length parameters
outer_min_length=150  # Minimum outer amplicon length
outer_max_length=500  # Maximum outer amplicon length
inner_min_length=80   # Minimum inner amplicon length
inner_max_length=300  # Maximum inner amplicon length

# Execute the semi-nested primer combinations script with the required file paths and amplicon length parameters as arguments
python "$semi_nested_script" "mapping_positions.tsv" "primer_metadata.tsv" $outer_min_length $outer_max_length $inner_min_length $inner_max_length

# Unload Python module
module unload python/3.9.15
