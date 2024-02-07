#!/bin/bash
#SBATCH --job-name=fetch_dengue_ncbi  # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH -p Priority,Background
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie,nelson_krisanto@hotmail.com # Where to send mail
#SBATCH --error="/home/people/23203786/scratch/Benchmarking_PCR_ONT/logs/fetch_dengue_ncbi_err_%j.log" # Standard error log ####
#SBATCH --output="/home/people/23203786/scratch/Benchmarking_PCR_ONT/logs/fetch_dengue_ncbi_%j.log" # Standard output log
#SBATCH --cpus-per-task=5

# Load configuration from a separate file
source "config.sh"

# Define the output directory for dengue virus sequences
output_dir="$my_dir/$project_dir/raw_data/NCBI_dengue_data" ## adjust
mkdir -p "$output_dir"

# Perform the search to get dengue virus sequences. Adjust the search term as needed.
esearch -db nuccore -query "Dengue Virus[ORGN] AND 300:11000[SLEN]" | efetch -format fasta > "$output_dir/dengue_virus_sequences.fasta"

echo "Fetched dengue virus sequences are saved in: $output_dir/dengue_virus_sequences.fasta"

# Specify the input FASTA file path for dengue virus
fasta_file="$output_dir/dengue_virus_sequences.fasta"

# Specify the output TSV file path for dengue virus sequence information
tsv_file="$output_dir/dengue_sequence_info.tsv"

# Rename fasta headers with accession IDs only
sed -i 's/^>\([^ ]*\) .*/>\1/' "$fasta_file"

# Extract accession IDs and sequence lengths to a TSV file with column headers
echo -e "accession_id\tsequence_length" > "$tsv_file"
awk -F '[>\t]' '/^>/ {if (acc != "") print acc "\t" len; acc=$2; len=0; next} {len+=length($0)} END {if (acc != "") print acc "\t" len}' "$fasta_file" >> "$tsv_file"

echo "Dengue virus sequence information is stored in: $tsv_file"

# Count the number of dengue virus sequences
num_sequences=$(grep -c '^>' "$fasta_file")
echo "Number of dengue virus sequences: $num_sequences"
