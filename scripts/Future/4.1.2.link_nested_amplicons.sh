#!/bin/bash
#SBATCH --job-name=4.1.2.link_nested_amplicons
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10

# Load Python module
module load python/3.9.15

# Define input file paths
outer_results_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/nested_outer_insilico_PCR_results.txt"
inner_results_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/nested_inner_insilico_PCR_results.txt"
inner_amplicons_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/nested_inner_insilico_PCR_amplicons.fasta"
primers_file="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/inner_nested_primer_combinations.tsv"
output_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/filtered_linked_nested_amplicons.csv"

# Path to the link_amplicons script
link_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/link_nested_amplicons.py"

# Run the linking script to filter amplicons between 100-500bp and link them to their accession numbers and primer pairs
python "$link_script" "$outer_results_file" "$inner_results_file" "$inner_amplicons_file" "$primers_file" "$output_file"

# Unload Python module
module unload python/3.9.15
