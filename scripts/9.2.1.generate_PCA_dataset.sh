#!/bin/bash
#SBATCH --job-name=9.1.1.generate_FAMD_dataset
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10

# Directory where the TSV files are located
tsv_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv"
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/conventional"

# Load Python module
module load python/3.9.15

# Change to the directory where the TSV files are located
cd "$tsv_dir"

# Path to the primer combinations script
combine_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/PCA_dataset.py"

# Define input and output file paths
mapping_file="$tsv_dir/mapping_positions.tsv"
metadata_file="$tsv_dir/primer_metadata.tsv"
combinations_file="$tsv_dir/conventional_primer_combinations.tsv"
output_file="$output_dir/individual_conv_primer_details.tsv"

# Execute the primer combinations script with the input and output files as arguments
python "$combine_script" "$mapping_file" "$metadata_file" "$combinations_file" "$output_file"

# Unload Python module
module unload python/3.9.15
    