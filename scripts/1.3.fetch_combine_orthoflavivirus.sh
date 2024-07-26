#!/bin/bash
#SBATCH --job-name=1.3.fetch_combine_orthoflavivirus  # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie,nelson_krisanto@hotmail.com # Where to send mail
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt  # Error log
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt  # Standard output log
#SBATCH --cpus-per-task=5

# Define paths
my_dir="/home/people/23203786/scratch/Nelson-Dissertation"
output_dir="$my_dir/raw_data/Orthoflavivirus"
logs_path="$my_dir/logs"
entrezdirect="/home/people/23203786/tools/entrez-direct_env"

# Load anaconda module
module load anaconda/3.5.2

# Activate your conda environment for Entrez Direct
source activate $entrezdirect

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Ensure the output directory was created successfully
if [ ! -d "$output_dir" ]; then
    echo "Failed to create output directory $output_dir. Exiting."
    exit 1
fi

# Define a list of accession numbers to download
accessions=(
    "DQ235145" "AY323490" "AF331718" "AF253419" "Y07863"
    "D12937" "X86784" "DQ235152" "DQ235151" "DQ235153"
    "AY193805" "L06436" "AF311056" "DQ235149" "U27495"
    "X07755" "L40361" "DQ235144" "DQ235150" "KF815939"
    "DQ235146" "AY632536" "AF013366" "AF013375" "AF013390"
    "U88536" "U87411" "M93130" "AF326573" "KF917536"
    "M18370" "AF013384" "AF161266" "DQ525916" "AY453411"
    "D00246" "M12294" "AF013413" "AY632541" "AF013407"
    "AY632545" "AY632539" "AF013397" "AF013377" "JX236040"
    "JF895923" "AY632535" "DQ837642" "EU707555" "X03700"
    "AY632540" "DQ859056" "DQ859057" "DQ859060" "DQ859066"
    "DQ859067" "DQ859062" "DQ859065" "DQ837641" "AF013405"
    "AB114858" "AF160193" "AF013370" "KJ469371" "AJ242984"
    "AF013401" "AF013402" "AF013365" "AF013368" "AF013371"
    "AJ299445" "AF013369" "AF013394" "AF144692"
)

# Convert the list to a comma-separated string
accession_list=$(IFS=, ; echo "${accessions[*]}")

# Fetch sequences from NCBI using the list of accession numbers
esearch -db nuccore -query "$accession_list" | efetch -format fasta > "$output_dir/orthoflavivirus.fasta"

echo "Fetched sequences saved in: $output_dir/orthoflavivirus.fasta"

# Specify the input FASTA file path for orthoflavivirus
fasta_file="$output_dir/orthoflavivirus.fasta"

# Specify the output TSV file path for orthoflavivirus sequence information
tsv_file="$output_dir/orthoflavivirus.tsv"

# Rename fasta headers with accession IDs only
sed -i 's/^>\([^ ]*\) .*/>\1/' "$fasta_file"

# Extract accession IDs and sequence lengths to a TSV file with column headers
echo -e "accession_id\tsequence_length" > "$tsv_file"
awk -F '[>\t]' '/^>/ {if (acc != "") print acc "\t" len; acc=$2; len=0; next} {len+=length($0)} END {if (acc != "") print acc "\t" len}' "$fasta_file" >> "$tsv_file"

echo "Orthoflavivirus sequence information is stored in: $tsv_file"

# Count the number of orthoflavivirus sequences
num_sequences=$(grep -c '^>' "$fasta_file")
echo "Number of orthoflavivirus sequences: $num_sequences"

# Deactivate the conda environment
conda deactivate


# Define paths
my_dir="/home/people/23203786/scratch/Nelson-Dissertation"
orthoflavivirus_fasta="$my_dir/raw_data/Orthoflavivirus/orthoflavivirus.fasta"
clustered_fasta="$my_dir/results/fasta/clustered_sequences.fasta"
combined_fasta="$my_dir/results/fasta/combined_sequences.fasta"

# Combine the two FASTA files
cat "$orthoflavivirus_fasta" "$clustered_fasta" > "$combined_fasta"

echo "Combined FASTA file created at: $combined_fasta"
