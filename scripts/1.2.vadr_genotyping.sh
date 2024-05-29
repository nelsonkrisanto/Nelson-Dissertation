#!/bin/bash
#SBATCH --job-name=1.3.vadr_genotyping
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Load the Perl module
module load perl/5.30 

# Define input and output file paths
input_fasta="/home/people/23203786/scratch/Nelson-Dissertation/results/mmseqs2_clustering/clustered_sequences.fasta"
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/vadr_genotyping"
mkdir -p "$output_dir"

# Run VADR
vadr_dir="/home/people/23203786/tools/vadr"
vadr_model_dir="$vadr_dir/vadr-models"
output_prefix="$output_dir/vadr_output"

perl $vadr_dir/v-annotate.pl --cpu 10 --mdir $vadr_model_dir --glsearch $input_fasta $output_prefix

# Unload the Perl module
module unload perl/5.30
