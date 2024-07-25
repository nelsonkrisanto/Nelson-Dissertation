#!/bin/bash
#SBATCH --job-name=2.1.fetch_ref_genome  # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie,nelson_krisanto@hotmail.com # Where to send mail
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt  # Error log
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt  # Standard output log
#SBATCH --cpus-per-task=2

# Load Anaconda module
module load anaconda/3.5.2

# Source conda.sh directly to initialize conda
source /opt/software/anaconda/3.5.2/etc/profile.d/conda.sh

# Activate conda environment with Entrez Direct installed
conda activate /home/people/23203786/tools/entrez-direct_env

# Define output directory and make sure it exists
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/references/tes"
mkdir -p "$output_dir"

# Fetch specific Dengue Virus reference genomes using their accession numbers
accessions=("GCF_000871845.1" "GCF_000862125.1" "GCF_000866625.1" "GCF_000865065.1")

# Create or clear the output fasta file
output_fasta="$output_dir/dengue_virus_ref_genomes.fasta"
> "$output_fasta"

for acc in "${accessions[@]}"
do
    echo "Fetching genome for accession: $acc"
    esearch -db assembly -query "$acc" | elink -target nucleotide | efetch -format fasta >> "$output_fasta"
    if [ $? -ne 0 ]; then
        echo "Failed to fetch genome for accession: $acc" >> "$output_dir/error_log.txt"
        exit 1
    fi
done

echo "Fetched Dengue Virus reference genomes saved to $output_fasta"

# Deactivate conda environment
conda deactivate
