#!/bin/bash
#SBATCH --job-name=02.fetch_papers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=2

# Set the path to your reference genomes
reference_path="/home/people/23203786/scratch/Nelson-Dissertation/references"

# Load the BWA module
module load bwa/0.7

# Iterate over each reference genome file and index it
for ref_file in "$reference_path"/*.fna; do
    echo "Indexing $ref_file..."
    bwa index "$ref_file"
    
    if [ $? -ne 0 ]; then
        echo "Failed to index $ref_file"
        # You can choose to exit the script or continue with the next file
        # exit 1
    else
        echo "Successfully indexed $ref_file"
    fi
done

echo "Indexing complete."
