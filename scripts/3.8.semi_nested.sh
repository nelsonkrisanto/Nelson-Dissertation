#!/bin/bash
#SBATCH --job-name=3.8.semi_nested
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10

# Directory where the TSV files are located
tsv_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv"

# Load Python module
module load python/3.9.15

# Change to the directory where the TSV files are located
cd "$tsv_dir"

# Path to the semi-nested primer combinations script (update this path to your actual semi-nested PCR script location)
semi_nested_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/semi_nested.py"

# Execute the semi-nested primer combinations script with the required file paths as arguments
python "$semi_nested_script" "mapping_positions.tsv" "primer_metadata.tsv"

# Unload Python module
module unload python/3.9.15
