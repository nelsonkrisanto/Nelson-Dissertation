#!/bin/bash
#SBATCH --job-name=1.2.cluster_MMseqs2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10

# Load the necessary modules (assuming MMseqs2 is available as a module)
module load mmseqs2

# Define input and output paths
INPUT_FASTA="/home/people/23203786/scratch/Nelson-Dissertation/raw_data/NCBI_dengue_data/1.1.dengue_virus_sequences.fasta"
OUTPUT_DIR="/home/people/23203786/scratch/Nelson-Dissertation/results/mmseqs2_clustering"
DB_DIR="$OUTPUT_DIR/db"
TMP_DIR="$OUTPUT_DIR/tmp"
CLUSTER_RESULT="$OUTPUT_DIR/clusters"

# Create output and temporary directories
mkdir -p $OUTPUT_DIR
mkdir -p $DB_DIR
mkdir -p $TMP_DIR

# Create MMseqs2 database from input fasta file
mmseqs createdb $INPUT_FASTA $DB_DIR/input_db

# Run MMseqs2 clustering with 97% identity threshold
mmseqs cluster $DB_DIR/input_db $CLUSTER_RESULT $TMP_DIR --min-seq-id 0.97 -c 0.8 --cov-mode 1 --threads 10

# Create output fasta file containing representative sequences of each cluster
mmseqs createsubdb $CLUSTER_RESULT $DB_DIR/input_db $OUTPUT_DIR/clustered_db
mmseqs convert2fasta $OUTPUT_DIR/clustered_db $OUTPUT_DIR/clustered_sequences.fasta

echo "MMseqs2 clustering completed. Results are in the $OUTPUT_DIR directory."

# Unload the MMseqs2 module
module unload mmseqs2
