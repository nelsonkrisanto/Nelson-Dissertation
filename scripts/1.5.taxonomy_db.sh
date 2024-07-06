#!/bin/bash
#SBATCH --job-name=1.5.taxonomy_db
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00

# Load Python module
module load python/3.7.4

# Define paths
#python_script_taxonomy="/home/people/23203786/scratch/Nelson-Dissertation/scripts/protein_db_tax.py"
combined_fasta_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/combine_fasta.py"

# Run the Python script to modify the taxonomy file
#python $python_script_taxonomy

echo "Python script completed. Taxonomy file saved to /home/people/23203786/scratch/Nelson-Dissertation/taxonomy/protein_db_taxonomy.tsv"

# Run the Python script to combine the FASTA files and rename headers
python $combined_fasta_script

echo "Combined FASTA files created."