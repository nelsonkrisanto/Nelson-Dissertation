#!/bin/bash
#SBATCH --job-name=3.6.insilico_PCR
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Load the Perl module
module load perl/5.30 

# Define input file paths
primers="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/conventional/top_200_conv_primer_combinations.tsv"
sequences="/home/people/23203786/scratch/Nelson-Dissertation/raw_data/NCBI_dengue_data/1.1.dengue_virus_sequences.fasta"
out_dir="/home/people/23203786/scratch/Nelson-Dissertation/results"
tool="/home/people/23203786/tools/in_silico_PCR/in_silico_PCR.pl" 

# Create output directory if it doesn't exist
mkdir -p "$out_dir"

# Change to the output directory
cd "$out_dir" || exit

# Run in-silico PCR
perl "$tool" -s "$sequences" -p "$primers" -m -i -r > insilico_PCR_results.txt 2> insilico_PCR_amplicons.fasta

# Unload the Perl module
module unload perl/5.30 
