#!/bin/bash
#SBATCH --job-name=blastn_classifier
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00


#!/bin/bash

# Create necessary directories
mkdir -p /home/people/23203786/scratch/Nelson-Dissertation/classified_blast
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

# Create directory for classified BLAST output
mkdir -p /home/people/23203786/scratch/Nelson-Dissertation/classified_blast

# Load Python module
module load python/3.7.4

# Clone the Assign-Taxonomy-with-BLAST repository if it doesn't exist
if [ ! -d "/home/people/23203786/scratch/Nelson-Dissertation/Assign-Taxonomy-with-BLAST" ]; then
    git clone https://github.com/Joseph7e/Assign-Taxonomy-with-BLAST.git /home/people/23203786/scratch/Nelson-Dissertation/Assign-Taxonomy-with-BLAST
fi

# Path to taxonomy assignment script
assign_taxonomy="/home/people/23203786/scratch/Nelson-Dissertation/Assign-Taxonomy-with-BLAST/taxonomy_assignment_BLAST.py"

# Ensure Biopython is installed
pip install biopython --user

# Run BLAST searches and assign taxonomy for each query FASTA file in the vsearch directory
for i in /home/people/23203786/scratch/Nelson-Dissertation/vsearch/*.fasta; do
    name=$(basename "$i" .fasta)

    # Run BLAST searches against each serotype database
    for serotype in denv1 denv2 denv3 denv4; do
        blastn -db /home/people/23203786/scratch/Nelson-Dissertation/blast_db/$serotype -query $i -out /home/people/23203786/scratch/Nelson-Dissertation/classified_blast/${serotype}_${name}.blastout.txt -num_threads 5 -outfmt "6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore"
        
        python3 $assign_taxonomy --blast_flavor blastn --cutoff_species 90 --cutoff_family 70 --cutoff_phylum 60 --length_percentage 0.50 --blast_database $serotype --output_dir /home/people/23203786/scratch/Nelson-Dissertation/classified_blast/${serotype}_${name} --blast_file /home/people/23203786/scratch/Nelson-Dissertation/classified_blast/${serotype}_${name}.blastout.txt $i /home/people/23203786/scratch/Nelson-Dissertation/taxonomy/custom_taxonomy.tsv
    done
done

# Extract and merge taxonomy results
find "/home/people/23203786/scratch/Nelson-Dissertation/classified_blast/" -type f -iname 'taxonomy_assignment_per_sequence.tsv' -exec sh -c '
    path="${1%/*}"; filename="${1##*/}";
    cp -nv "${1}" "/home/people/23203786/scratch/Nelson-Dissertation/classified_blast/${path##*/}.tsv" ' sh_cp {} \;

# Remove subdirectories and content
rm -rf /home/people/23203786/scratch/Nelson-Dissertation/classified_blast/*/

# Merge all taxonomy outputs but include file name
awk 'FNR>1 || NR==1 {print $0","FILENAME}' /home/people/23203786/scratch/Nelson-Dissertation/classified_blast/*.tsv > /home/people/23203786/scratch/Nelson-Dissertation/classified_blast/taxonomy_blast.tsv

echo "BLAST search and taxonomy assignment completed successfully."
