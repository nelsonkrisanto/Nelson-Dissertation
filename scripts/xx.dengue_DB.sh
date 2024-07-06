#!/bin/bash
#SBATCH --job-name=1.3.diamond_blastx  # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie,nelson_krisanto@hotmail.com # Where to send mail
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_diamond_%x_%j.txt  # Error log
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_diamond_%x_%j.txt  # Standard output log
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Load configuration from a separate file
my_dir="/home/people/23203786/scratch/Nelson-Dissertation"

# Load anaconda module for DIAMOND
module load anaconda/3.5.2

# Activate your conda environment for DIAMOND
echo "Activating DIAMOND environment..."
source activate /home/people/23203786/tools/diamond_env

# Verify environment activation
if [ $? -ne 0 ]; then
    echo "Failed to activate DIAMOND environment"
    exit 1
fi

# Define paths
fasta_file="$my_dir/raw_data/tes/dengue_virus_sequences.fasta"
db_output="/home/people/23203786/scratch/Nelson-Dissertation/results/diamond_db/test_dengue_db"
tmp_dir="/home/people/23203786/scratch/Nelson-Dissertation/tmp"
output_file="${db_output}/matches.m8"

# Create necessary directories
mkdir -p "$tmp_dir"
mkdir -p "$(dirname "$db_output")"
mkdir -p "$db_output"

# Verify directory creation
echo "Temporary directory: $tmp_dir"
echo "Database output directory: $db_output"
echo "Output file path: $output_file"

# Build DIAMOND database
echo "Building DIAMOND database..."
diamond makedb --in "$fasta_file" -d "$db_output"

# Check if the database was built successfully
if [ $? -ne 0 ]; then
    echo "Failed to build DIAMOND database"
    exit 1
fi

# Run DIAMOND blastx
echo "Running DIAMOND blastx..."
diamond blastx -d "$db_output" -q "$fasta_file" -o "$output_file" --outfmt 6 --evalue 1 --max-target-seqs 1000 --tmpdir "$tmp_dir" --verbose

# Check for the DIAMOND blastx output
if [ ! -f "$output_file" ]; then
    echo "Failed to run DIAMOND blastx"
    exit 1
fi

echo "DIAMOND blastx completed successfully."
cat "$output_file"

# Deactivate Conda environment
conda deactivate
