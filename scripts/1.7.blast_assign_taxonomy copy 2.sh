#!/bin/bash
#SBATCH --job-name=1.7.blast_and_assign_taxonomy
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Load necessary modules
module load blast/2.12.0
source activate /home/people/23203786/tools/biopython

# Define paths
query_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/passed_sequences.fasta"
unique_query_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/unique_passed_sequences.fasta"
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/blast_output"
taxonomy_script="/home/people/23203786/scratch/Assign-Taxonomy-with-BLAST/taxonomy_assignment_BLAST.py"
taxonomy_file="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/protein_db_taxonomy.tsv"

# Ensure the output directory exists and is empty
if [ -d "$output_dir" ]; then
    echo "Output directory exists. Clearing old files."
    rm -r $output_dir/*
else
    mkdir -p $output_dir
fi

# Download and extract seqkit if not already done
seqkit_path="/home/people/23203786/tools/seqkit/seqkit"
if [ ! -f "$seqkit_path" ]; then
    echo "Downloading seqkit..."
    wget -O /home/people/23203786/tools/seqkit_linux_amd64.tar.gz https://github.com/shenwei356/seqkit/releases/download/v0.16.1/seqkit_linux_amd64.tar.gz
    tar -xzvf /home/people/23203786/tools/seqkit_linux_amd64.tar.gz -C /home/people/23203786/tools/
    seqkit_path="/home/people/23203786/tools/seqkit"
fi

# Remove duplicate sequences
echo "Removing duplicate sequences..."
$seqkit_path/seqkit rmdup -s -i $query_file -o $unique_query_file

# Define serotypes and their database paths
declare -A serotypes=(
    ["denv4"]="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/denv4/dengue_blast_db"
)

for serotype in "${!serotypes[@]}"; do
    db_path="${serotypes[$serotype]}"
    blast_output="${output_dir}/${serotype}_blast_output.txt"
    
    # Run BLAST search
    echo "Running BLAST search for $serotype..."
    blastx -db $db_path -query $unique_query_file -out $blast_output -evalue 1e-5 -outfmt "6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore" -num_threads 10
    
    # Assign taxonomy
    taxonomy_output="${output_dir}/${serotype}_taxonomy_assignment.tsv"
    echo "Assigning taxonomy for $serotype..."
    python3 $taxonomy_script --blast_flavor blastx --cutoff_species 90 --cutoff_family 70 --cutoff_phylum 60 --length_percentage 0.50 --blast_database $db_path --output_dir $output_dir --blast_file $blast_output $unique_query_file $taxonomy_file
    
    echo "BLAST search and taxonomy assignment for $serotype completed."
done

# Unload modules
module unload blast/2.12.0
source deactivate

echo "All BLAST searches and taxonomy assignments completed successfully."
