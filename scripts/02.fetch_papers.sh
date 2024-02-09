#!/bin/bash
#SBATCH --job-name=02.fetch_papers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie,nelson_krisanto@hotmail.com # Where to send mail
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt  # Error log
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt  # Standard output log
#SBATCH --cpus-per-task=5

# Load configuration from a separate file
source "00.config.sh"

# Define variables
output_directory="$my_dir/data/primers"
query="https://pubmed.ncbi.nlm.nih.gov/?term=dengue+virus+primers"

# Create output directory
mkdir -p "$output_directory" || { echo "Error: Unable to create output directory."; exit 1; }

# Run Pubmed query using Python script
python $scrappaper --url "https://pubmed.ncbi.nlm.nih.gov/?term=dengue+virus+primers" --pages 40 > $output_directory/dengue_virus.tsv
