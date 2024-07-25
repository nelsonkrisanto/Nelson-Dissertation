#!/bin/bash
#SBATCH --job-name=3.6.2.nested_insilico_PCR
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=96:00:00

# Load the Perl module
module load perl/5.30 

# Define input file paths
outer_primers="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/outer_nested_primer_combinations.tsv"
nested_primers="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/inner_nested_primer_combinations.tsv"
sequences="/home/people/23203786/scratch/Nelson-Dissertation/results/fasta/combined_sequences.fasta"
out_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico"
tool="/home/people/23203786/tools/in_silico_PCR/in_silico_PCR.pl"
filter_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/filter_nested_amplicons.py"

# Log the start of the script
echo "Script started at $(date)"

# Change working directory
echo "Changing directory to $out_dir..."
cd "$out_dir" || exit

# Run in-silico PCR with outer primers
echo "Running in-silico PCR with outer primers..."
perl "$tool" -s "$sequences" -p "$outer_primers" -l 2000 -m -i -e -r > nested_outer_insilico_PCR_results.txt 2> nested_outer_insilico_PCR_amplicons.fasta
echo "Outer in-silico PCR completed."

# Run in-silico PCR with nested primers on the outer PCR products
echo "Running in-silico PCR with nested primers..."
perl "$tool" -s "nested_outer_insilico_PCR_amplicons.fasta" -p "$nested_primers" -l 500 -m -i -e -r > nested_inner_insilico_PCR_results.txt 2> nested_inner_insilico_PCR_amplicons.fasta
echo "Nested in-silico PCR completed."

# Unload Perl module
echo "Unloading Perl module..."
module unload perl/5.30

# Load Python module
echo "Loading Python module..."
module load python/3.9.15

# Run the filtering script to filter inner nested amplicons between 100-500bp
echo "Running filtering script..."
python "$filter_script" nested_inner_insilico_PCR_results.txt nested_inner_insilico_PCR_amplicons.fasta filtered_nested_insilico_PCR_amplicons.fasta
echo "Filtering script completed."

# Unload Python module
echo "Unloading Python module..."
module unload python/3.9.15

# Log the end of the script
echo "Script completed at $(date)"
