#!/bin/bash
#SBATCH --job-name=3.7.insilico_PCR
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10

module load perl/5.30 
# Load configuration from a separate file
source "$my_dir/$project_dir/config.sh"

# Input file paths
primers="$output_dir/insilico_pcr_primers.tsv"
sequences="$output_dir/test_norovirus.fasta"
out_dir="$my_dir/$project_dir/insilico_PCR"
tool="/home/people/fitzpatria/tools/in_silico_PCR.pl" 

# Create output directory
mkdir -p "$out_dir"

# Change working directory
cd "$out_dir" || exit

# Run in-silico PCR for VP1
perl "$tool" l -s "$sequences" -p "$primers" -m -i -r > insilico_PCR_results1.txt 2> insilico_PCR_amplicons1.fasta

module unload perl/5.30 