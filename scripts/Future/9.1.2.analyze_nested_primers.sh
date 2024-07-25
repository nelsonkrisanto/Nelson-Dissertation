#!/bin/bash
#SBATCH --job-name=9.1.2.analyze_nested_primers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10

# Directory where the TSV files are located
tsv_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv"
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv"

# Load Python module
module load python/3.9.15

# Change to the directory where the TSV files are located
cd "$tsv_dir"

# Path to the primer combinations script
combine_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/analyze_primers.py"

# Define input and output file paths
input_file="$tsv_dir/mapping_positions.tsv"
output_file="$output_dir/nested/individual_nested_primer_details.tsv"

# Execute the primer combinations script with the input and output files as arguments
python "$combine_script" "$input_file" "$output_file"

# Unload Python module
module unload python/3.9.15
