#!/bin/bash
#SBATCH --job-name=3.5.3.semi_nested
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

# Path to the nested primer combinations script
nested_pcr_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/semi_nested.py"

# Execute the nested primer combinations script with the necessary arguments
python "$nested_pcr_script" "mapping_positions.tsv" "primer_metadata.tsv" 300 100

# Unload Python module
module unload python/3.9.15
