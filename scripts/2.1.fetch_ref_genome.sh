#!/bin/bash
#SBATCH --job-name=2.1.fetch_ref_genome  # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie,nelson_krisanto@hotmail.com # Where to send mail
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt  # Error log
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt  # Standard output log
#SBATCH --cpus-per-task=2

# Load modules (if required on your cluster)
module load anaconda/3.5.2

# Activate conda environment with Entrez Direct installed
source activate entrez_direct_env

# Define output directory and make sure it exists
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/references/tes"
mkdir -p "$output_dir"

# Fetch Dengue Virus reference genomes
# Note: Adjust the query as needed to refine the results
esearch -db nucleotide -query "Dengue virus[Organism] AND refseq[Filter]" | \
efetch -format fasta > "$output_dir/dengue_virus_ref_genomes.fasta"

echo "Fetched Dengue Virus reference genomes saved to $output_dir/dengue_virus_ref_genomes.fasta"

# Deactivate conda environment
conda deactivate
