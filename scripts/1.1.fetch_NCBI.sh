#!/bin/bash
#SBATCH --job-name=1.1.fetch_NCBI  # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie,nelson_krisanto@hotmail.com # Where to send mail
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt  # Error log
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt  # Standard output log
#SBATCH --cpus-per-task=5

# Load configuration from a separate file
source "00.config.sh"

# Load anaconda module
module load anaconda/3.5.2

# Activate your conda environment for Entrez Direct
source activate $entrezdirect

# Define the output directory for dengue virus sequences
output_dir="$my_dir/raw_data/NCBI_dengue_data" ## adjust
mkdir -p "$output_dir"

# Perform the search to get dengue virus sequences. Adjust the search term as needed.
esearch -db nuccore -query "Dengue Virus[ORGN] AND 300:11000[SLEN]" | efetch -format fasta > "$output_dir/dengue_virus_sequences.fasta" ## adjust

echo "Fetched dengue virus sequences are saved in: $output_dir/dengue_virus_sequences.fasta" ## adjust

# Specify the input FASTA file path for dengue virus
fasta_file="$output_dir/dengue_virus_sequences.fasta" ## adjust

# Specify the output TSV file path for dengue virus sequence information
tsv_file="$output_dir/dengue_sequence_info.tsv" ## adjust

# Rename fasta headers with accession IDs only
sed -i 's/^>\([^ ]*\) .*/>\1/' "$fasta_file"

# Extract accession IDs and sequence lengths to a TSV file with column headers
echo -e "accession_id\tsequence_length" > "$tsv_file"
awk -F '[>\t]' '/^>/ {if (acc != "") print acc "\t" len; acc=$2; len=0; next} {len+=length($0)} END {if (acc != "") print acc "\t" len}' "$fasta_file" >> "$tsv_file"

echo "Dengue virus sequence information is stored in: $tsv_file" ## adjust

# Count the number of dengue virus sequences
num_sequences=$(grep -c '^>' "$fasta_file")
echo "Number of dengue virus sequences: $num_sequences"

# Deactivate the conda environment
conda deactivate