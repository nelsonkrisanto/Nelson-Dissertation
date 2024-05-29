#!/bin/bash
#SBATCH --job-name=3.6.2.nested_insilico_PCR
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Load the Perl module
module load perl/5.30 

# Define input file paths
outer_primers="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/nested_outer_primer_combinations.tsv"
inner_primers="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/nested_inner_primer_combinations.tsv"
sequences="/home/people/23203786/scratch/Nelson-Dissertation/raw_data/NCBI_dengue_data/1.1.dengue_virus_sequences.fasta"
out_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico"
tool="/home/people/23203786/tools/in_silico_PCR/in_silico_PCR.pl"

# Create output directory if it doesn't exist
mkdir -p "$out_dir"

# Change working directory
cd "$out_dir" || exit

# Make sure the Perl script is executable
chmod +x "$tool"

# Run in-silico PCR for outer primers
perl "$tool" -s "$sequences" -p "$outer_primers" -m 2 -i -e -r > nested_outer_insilico_PCR_results.txt 2> nested_outer_insilico_PCR_amplicons.fasta

# Check if the outer PCR results were successfully generated
if [ $? -eq 0 ]; then
  # Use the outer amplicons as input for the second round
  perl "$tool" -s nested_outer_insilico_PCR_amplicons.fasta -p "$inner_primers" -m 2 -i -e -r > nested_inner_insilico_PCR_results.txt 2> nested_inner_insilico_PCR_amplicons.fasta

  # Check if the inner PCR results were successfully generated
  if [ $? -ne 0 ]; then
    echo "Error: Inner in-silico PCR failed."
    exit 1
  fi
else
  echo "Error: Outer in-silico PCR failed."
  exit 1
fi

# Unload the Perl module
module unload perl/5.30
