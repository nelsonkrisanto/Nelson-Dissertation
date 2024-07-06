#!/bin/bash
#SBATCH --job-name=1.6.convert_and_create_blast_db
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Load necessary modules
module load python/3.7.4

# Set PYTHONPATH to include FROGS library
export PYTHONPATH=/home/people/23203786/scratch/FROGS/lib:$PYTHONPATH

SEQKIT="/home/people/23203786/tools/seqkit"

# Check available BLAST versions and use an alternative if necessary
module avail blast

# Set the appropriate BLAST module version (update this based on the output of the module avail command)
module load blast/2.12.0

# Define paths
python_script="/home/people/23203786/scratch/FROGS/libexec/fasta2RDP.py"
input_taxonomy="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/protein_db_taxonomy.tsv"
input_fasta="/home/people/23203786/scratch/Nelson-Dissertation/db/denv4/combined_denv4.fasta"
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/denv4"

# Create output directory if it doesn't exist
mkdir -p $output_dir
output_rdp_fasta="$output_dir/protein_db_taxonomy_RDP.fasta"
output_rdp_taxonomy="$output_dir/protein_db_taxonomy_RDP.tax"

# Remove spaces from FASTA headers using seqkit
echo "Removing spaces from FASTA headers..."
$SEQKIT replace -p " " -r "_" -o $output_dir/cleaned_denv4.fasta $input_fasta

# Run the Python script to convert to RDP format
echo "Converting denv4 to RDP format..."
python3.7 $python_script -d $output_dir/cleaned_denv4.fasta -t $input_taxonomy -r Domain Phylum Class Order Family Genus Species --rdp-taxonomy $output_rdp_taxonomy --rdp-fasta $output_rdp_fasta

# Check if the conversion was successful
if [[ $? -ne 0 ]]; then
    echo "Conversion to RDP format failed for denv4"
    exit 1
fi

# Validate the converted FASTA file for invalid residues
echo "Validating FASTA file for denv4..."
grep -E -v '^[A-Z]+$|^>' $output_rdp_fasta
if [[ $? -eq 0 ]]; then
    echo "FASTA file for denv4 contains invalid residues."
    exit 1
fi

# Create BLAST database
echo "Creating BLAST database for denv4..."
makeblastdb -in $output_rdp_fasta -parse_seqids -taxid_map $output_rdp_taxonomy -dbtype prot -out $output_dir/dengue_blast_db

# Check if the BLAST database creation was successful
if [[ $? -ne 0 ]]; then
    echo "BLAST database creation failed for denv4"
    exit 1
fi

echo "BLAST database for denv4 created successfully."

# Unload modules
module unload python/3.7.4
module unload blast/2.12.0

echo "All BLAST databases created successfully."
