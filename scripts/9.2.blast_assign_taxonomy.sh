#!/bin/bash
#SBATCH --job-name=blastassigntax
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Create necessary directories
mkdir -p /home/people/23203786/scratch/Nelson-Dissertation/classified_blast

# Load BLAST module
module load blast/2.12.0

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

# Verify query FASTA files format in vsearch directory
for query_file in /home/people/23203786/scratch/Nelson-Dissertation/vsearch/*.fasta; do
    echo "Verifying format of $query_file"
    head -n 2 $query_file
done

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
