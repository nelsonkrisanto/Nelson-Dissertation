#!/bin/bash
#SBATCH --job-name=blastdb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Create necessary directories
mkdir -p /home/people/23203786/scratch/Nelson-Dissertation/blast_db

# Load BLAST module
module load blast/2.12.0

# Copy necessary files
cp /home/people/23203786/scratch/Nelson-Dissertation/db/denv1/denv1.fasta /home/people/23203786/scratch/Nelson-Dissertation/blast_db/denv1.fasta
cp /home/people/23203786/scratch/Nelson-Dissertation/db/denv2/denv2.fasta /home/people/23203786/scratch/Nelson-Dissertation/blast_db/denv2.fasta
cp /home/people/23203786/scratch/Nelson-Dissertation/db/denv3/denv3.fasta /home/people/23203786/scratch/Nelson-Dissertation/blast_db/denv3.fasta
cp /home/people/23203786/scratch/Nelson-Dissertation/db/denv4/denv4.fasta /home/people/23203786/scratch/Nelson-Dissertation/blast_db/denv4.fasta

cp /home/people/23203786/scratch/Nelson-Dissertation/taxonomy/denv1.tax /home/people/23203786/scratch/Nelson-Dissertation/blast_db/denv1.tax
cp /home/people/23203786/scratch/Nelson-Dissertation/taxonomy/denv2.tax /home/people/23203786/scratch/Nelson-Dissertation/blast_db/denv2.tax
cp /home/people/23203786/scratch/Nelson-Dissertation/taxonomy/denv3.tax /home/people/23203786/scratch/Nelson-Dissertation/blast_db/denv3.tax
cp /home/people/23203786/scratch/Nelson-Dissertation/taxonomy/denv4.tax /home/people/23203786/scratch/Nelson-Dissertation/blast_db/denv4.tax

cp /home/people/23203786/scratch/Nelson-Dissertation/taxonomy/custom_taxonomy.tsv /home/people/23203786/scratch/Nelson-Dissertation/blast_db/custom_taxonomy.tsv

cd /home/people/23203786/scratch/Nelson-Dissertation/blast_db

# Create BLAST databases for each input FASTA file
for i in *.fasta; do
    name=$(basename "$i" .fasta)
    makeblastdb -in $i -parse_seqids -taxid_map ${name}.tax -dbtype nucl -out ${name}
done

echo "BLAST databases created successfully."
