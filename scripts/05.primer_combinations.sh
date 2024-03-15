#!/bin/bash
#SBATCH --job-name=02.fetch_papers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=5

# Load configuration from a separate file
#source "$my_dir/$project_dir/config.sh"

# Load Python module
module load python/3.9.15

# Path to the primer combinations script
combine_script="$my_dir/$project_dir/scripts/primer_combinations.py"

# Execute the primer combinations script
python "$combine_script" "/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/mapping_positions.tsv" "/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/primer_metadata.tsv"

# Unload Python module
module unload python/3.9.15
