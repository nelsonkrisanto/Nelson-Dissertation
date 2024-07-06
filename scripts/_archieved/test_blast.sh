#!/bin/bash
#SBATCH --job-name=testblast
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Load BLAST module
module load blast/2.12.0

# Define paths
blast_db_dir="/home/people/23203786/scratch/Nelson-Dissertation/blast_db"
taxonomy_dir="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy"
mkdir -p $blast_db_dir

# Create BLAST databases
for serotype in denv1 denv2 denv3 denv4; do
    fasta_file="/home/people/23203786/scratch/Nelson-Dissertation/db/$serotype/$serotype.fasta"
    taxid_map_path="$taxonomy_dir/${serotype}.tax"
    makeblastdb -in $fasta_file -parse_seqids -taxid_map $taxid_map_path -dbtype nucl -out $blast_db_dir/$serotype
done

echo "BLAST databases created successfully."