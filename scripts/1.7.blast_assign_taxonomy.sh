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
base_output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/blast_output"
taxonomy_script="/home/people/23203786/scratch/Assign-Taxonomy-with-BLAST/taxonomy_assignment_BLAST.py"
taxonomy_file="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/protein_db_taxonomy.tsv"
seqkit_path="/home/people/23203786/scratch/Nelson-Dissertation/tools/seqkit"

# Ensure the base output directory exists
mkdir -p $base_output_dir

# Clear old temporary directories (Optional)
rm -rf ${base_output_dir}/tmp.*

# Generate a unique directory for this run
output_dir=$(mktemp -d -p $base_output_dir)

# Remove duplicate sequences
echo "Removing duplicate sequences..."
$seqkit_path rmdup -s -i $query_file -o $unique_query_file

# Define serotypes and their database paths
declare -A serotypes=(
    ["denv4"]="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/denv4/dengue_blast_db"
)

for serotype in "${!serotypes[@]}"; do
    db_path="${serotypes[$serotype]}"
    blast_output="${output_dir}/${serotype}_blast_output.txt"
    
    # Run BLAST search
    echo "Running BLAST search for $serotype..."
    blastx -db $db_path -query $unique_query_file -out $blast_output -evalue 1e-5 -outfmt "6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore staxids" -num_threads 10

    # Assign taxonomy
    taxonomy_output="${output_dir}/${serotype}_taxonomy_assignment.tsv"
    echo "Assigning taxonomy for $serotype..."

    # Ensure output directory exists
    if [ -d "$output_dir" ]; then
        echo "Output directory already exists. Clearing old files."
        rm -rf ${output_dir}/*
    else
        mkdir -p $output_dir
    fi

    python3 $taxonomy_script --blast_flavor blastx --cutoff_species 90 --cutoff_family 70 --cutoff_phylum 60 --length_percentage 0.50 --blast_database $db_path --output_dir $output_dir --blast_file $blast_output $unique_query_file $taxonomy_file
    
    if [ $? -ne 0 ]; then
        echo "Taxonomy assignment failed for $serotype."
        exit 1
    fi

    echo "BLAST search and taxonomy assignment for $serotype completed."
done

# Unload modules
module unload blast/2.12.0
source deactivate

echo "All BLAST searches and taxonomy assignments completed successfully."
