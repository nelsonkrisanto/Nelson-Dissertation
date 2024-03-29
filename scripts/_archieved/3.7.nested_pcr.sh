#!/bin/bash
#SBATCH --job-name=3.7.nested_pcr
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

# Path to the primer combinations script
nested_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/nested_pcr.py"

# Execute the primer combinations script with the required file paths as arguments
python "$nested_script" "mapping_positions.tsv" "primer_metadata.tsv"

# Unload Python module
module unload python/3.9.15
