#!/bin/bash
#SBATCH --job-name=1.7.3.blast_assign_taxonomy_amp_semi_nested
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
query_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/passed_semi_nested_sequences.fasta"
unique_query_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/unique_passed_semi_nested_sequences.fasta"  # Path for unique sequences
base_output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/blast_output/amp_semi_nested"
taxonomy_script="/home/people/23203786/scratch/Assign-Taxonomy-with-BLAST/taxonomy_assignment_BLAST.py"
taxonomy_file="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/protein_db_taxonomy.tsv"
best_taxonomy_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/best_taxonomy_assignment.py"
seqkit_path="/home/people/23203786/scratch/Nelson-Dissertation/tools/seqkit"

# Create a temporary directory and set it as TMPDIR
TMP_DIR=$(mktemp -d)
export TMPDIR=$TMP_DIR
trap 'rm -rf "$TMP_DIR"' EXIT

# Ensure the base output directory exists
mkdir -p $base_output_dir

# Remove duplicate sequences
echo "Removing duplicate sequences..."
$seqkit_path rmdup -s -i $query_file -o $unique_query_file

# Define serotypes and their database paths
declare -A serotypes=(
    ["denv4"]="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/denv4/dengue_blast_db"
    ["denv3"]="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/denv3/dengue_blast_db"
    ["denv2"]="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/denv2/dengue_blast_db"
    ["denv1"]="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/denv1/dengue_blast_db"
)

# Run BLAST and taxonomy assignment for each serotype
for serotype in "${!serotypes[@]}"; do
    db_path="${serotypes[$serotype]}"
    blast_output="${base_output_dir}/${serotype}_blast_output.txt"
    
    # Run BLAST search
    echo "Running BLAST search for $serotype..."
    blastx -db $db_path -query $unique_query_file -out $blast_output -evalue 1e-5 -outfmt "6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore staxids" -num_threads 10

    # Check if BLAST output file is generated
    if [ ! -f "$blast_output" ]; then
        echo "BLAST search for $serotype failed. No output file generated."
        exit 1
    else
        echo "BLAST search for $serotype completed successfully. Output file: $blast_output"
    fi

    # Run taxonomy assignment
    echo "Running taxonomy assignment for $serotype..."
    python3 $taxonomy_script --blast_flavor blastx --cutoff_species 70 --cutoff_family 60 --cutoff_phylum 50 --output_dir $base_output_dir --verbose --blast_file $blast_output -- $unique_query_file $taxonomy_file
    if [ $? -ne 0 ]; then
        echo "Taxonomy assignment for $serotype failed."
        exit 1
    else
        echo "Taxonomy assignment for $serotype completed successfully."
    fi
done

# Combine taxonomy assignments and determine the best one for each sequence
echo "Determining the best taxonomy assignment for each sequence..."
combined_taxonomy_file="${base_output_dir}/taxonomy_assignment_per_sequence.tsv"
best_output_file="${base_output_dir}/best_taxonomy_assignment_per_sequence.tsv"

python3 $best_taxonomy_script $combined_taxonomy_file $best_output_file

# Unload modules
module unload blast/2.12.0
source deactivate

echo "All BLAST searches and taxonomy assignments completed successfully."
echo "Best taxonomy assignment per sequence has been saved to ${best_output_file}"