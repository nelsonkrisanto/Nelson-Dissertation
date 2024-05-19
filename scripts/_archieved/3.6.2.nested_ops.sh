#!/bin/bash
#SBATCH --job-name=3.6.2.nested_ops
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
nested_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/nested_ops.py"

# Define outer and inner amplicon lengths
outer_min_len=150
outer_max_len=500
inner_min_len=80
inner_max_len=300

# Execute the primer combinations script with the required file paths and amplicon length arguments
python "$nested_script" "mapping_positions.tsv" "primer_metadata.tsv" $outer_min_len $outer_max_len $inner_min_len $inner_max_len

# Unload Python module
module unload python/3.9.15
