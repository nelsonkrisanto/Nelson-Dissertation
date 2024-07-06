#!/bin/bash
#SBATCH --job-name=1.7.make_blast_db
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00

# Load BLAST module
module load blast/2.8.1

# Create BLAST database
makeblastdb -in /home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/protein_db_taxonomy_RDP.fasta \
  -parse_seqids -taxid_map /home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/protein_db_taxonomy_RDP.tax \
  -dbtype nucl -out /home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/dengue_blast_db

# Unload BLAST module
module unload blast/2.8.1

echo "BLAST database creation completed."
