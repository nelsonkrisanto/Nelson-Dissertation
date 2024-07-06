#!/bin/bash
#SBATCH --job-name=1.6.blastn
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Define paths
blast_db="/home/people/23203786/scratch/Nelson-Dissertation/diamond_db/dengue_nucleotide_db"
input_fasta="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/passed_sequences.fasta"
output_file="/home/people/23203786/scratch/Nelson-Dissertation/results/blastn_output.txt"

# Load BLAST module
module load blast/2.11.0

# Run BLASTN
blastn -db $blast_db -query $input_fasta -out $output_file -evalue 1e-5 -outfmt 6 -num_threads 10

echo "BLASTN completed. Results saved to $output_file"

# Unload BLAST module
module unload blast/2.11.0
