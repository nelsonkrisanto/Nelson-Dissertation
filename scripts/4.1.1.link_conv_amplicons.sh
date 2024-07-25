#!/bin/bash
#SBATCH --job-name=4.1.1.link_conv_amplicons
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10

# Load Python module
module load python/3.9.15

# Define input file paths
results_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/conv_insilico_PCR_results.txt"
amplicons_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/conv_insilico_PCR_amplicons.fasta"
primers_file="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/conventional_primer_combinations.tsv"
output_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/filtered_linked_amplicons.csv"

# Path to the link_amplicons script
link_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/link_conv_amplicons.py"

# Run the linking script to filter amplicons between 300-1200bp and link them to their accession numbers and primer pairs
python "$link_script" "$results_file" "$amplicons_file" "$primers_file" "$output_file"

# Unload Python module
module unload python/3.9.15
