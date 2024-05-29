#!/bin/bash
#SBATCH --job-name=1.3.cd_hit_clustering
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10

# Load the necessary modules if required (example: module load cd-hit)
# module load cd-hit

# Define input and output paths
INPUT_FASTA="/home/people/23203786/scratch/Nelson-Dissertation/results/mmseqs2_clustering/clustered_sequences.fasta"
OUTPUT_FASTA="/home/people/23203786/scratch/Nelson-Dissertation/results/clustered_refseq.fasta"

# Run CD-HIT with 97% identity threshold
/home/people/23203786/tools/cdhit/cd-hit -i $INPUT_FASTA -o $OUTPUT_FASTA -c 0.97 -n 5 -T 10 -M 16000

echo "CD-HIT clustering completed. Results are in $OUTPUT_FASTA"
