#!/bin/bash
#SBATCH --job-name=3.7.2.nested_combine_primers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=100:00:00

# Load the Python module
module load python/3.9.15

# Define input file paths
outer_results_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/nested_outer_insilico_PCR_results.txt"
inner_results_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/nested_inner_insilico_PCR_results.txt"
combined_output_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/combined_nested_insilico_PCR_results.tsv"

# Path to the Python script
combine_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/nested_combine_primers.py"

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$combined_output_file")"

# Change working directory
cd "$(dirname "$combined_output_file")" || exit

# Run the Python script
python "$combine_script" "$outer_results_file" "$inner_results_file" "$combined_output_file"

# Unload the Python module
module unload python/3.9.15
