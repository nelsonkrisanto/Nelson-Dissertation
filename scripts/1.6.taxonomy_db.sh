#!/bin/bash
#SBATCH --job-name=1.6.taxonomy_db
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00

# Load Python module
module load python/3.7.4

# Define paths
python_script_taxonomy="/home/people/23203786/scratch/Nelson-Dissertation/scripts/protein_db_tax.py"
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy"
combined_dir="/home/people/23203786/scratch/Nelson-Dissertation/db"

# Run the Python script to modify the taxonomy file
python $python_script_taxonomy

# Dictionary with strain names and their corresponding new lineage nomenclature
declare -A dengue_strains
dengue_strains=(
    ["Dengue virus 1 Brazil/97-11/1997"]="DENV-1_genotype_V"
    ["Dengue virus 1 Jamaica/CV1636/1977"]="DENV-1_genotype_II"
    ["Dengue virus 1 Nauru/West Pac/1974"]="DENV-1_genotype_I"
    ["Dengue virus 1 Singapore/S275/1990"]="DENV-1_genotype_IV"
    ["Dengue virus 1 Thailand/AHF 82-80/1980"]="DENV-1_genotype_III"
    ["Dengue virus 2 16681-PDK53"]="DENV-2_genotype_Asian_I"
    ["Dengue virus 2 China/D2-04"]="DENV-2_genotype_Asian_II"
    ["Dengue virus 2 Jamaica/1409/1983"]="DENV-2_genotype_American"
    ["Dengue virus 2 Peru/IQT2913/1996"]="DENV-2_genotype_Cosmopolitan"
    ["Dengue virus 3 China/80-2/1980"]="DENV-3_genotype_III"
    ["Dengue virus 3 Martinique/1243/1999"]="DENV-3_genotype_IV"
    ["Dengue virus 3 Philippines/H87/1956"]="DENV-3_genotype_I"
    ["Dengue virus 3 Singapore/8120/1995"]="DENV-3_genotype_II"
    ["Dengue virus 3 Sri Lanka/1266/2000"]="DENV-3_genotype_V"
    ["Dengue virus 4 Dominica/814669/1981"]="DENV-4_genotype_II"
    ["Dengue virus 4 Philippines/H241/1956"]="DENV-4_genotype_I"
)

# Base URL for downloading sequences
base_url="https://www.ncbi.nlm.nih.gov/protein/"

# Create combined directories if they don't exist
mkdir -p $combined_dir/denv1
mkdir -p $combined_dir/denv2
mkdir -p $combined_dir/denv3
mkdir -p $combined_dir/denv4

# Download, rename, and combine sequences
for strain in "${!dengue_strains[@]}"; do
    lineage=${dengue_strains[$strain]}
    serotype=$(echo $lineage | cut -d'_' -f1)
    
    # Extract accession number from the protein name
    accession=$(echo $strain | grep -oP 'P\w+')
    
    # Download the sequence
    wget -O ${lineage}.fasta ${base_url}${accession}.fasta
    
    # Rename the sequence header
    sed -i "1s/.*/>${lineage}/" ${lineage}.fasta
    
    # Combine sequences into respective files
    if [[ $serotype == "DENV-1" ]]; then
        cat ${lineage}.fasta >> $combined_dir/denv1/combined_denv1.fasta
    elif [[ $serotype == "DENV-2" ]]; then
        cat ${lineage}.fasta >> $combined_dir/denv2/combined_denv2.fasta
    elif [[ $serotype == "DENV-3" ]]; then
        cat ${lineage}.fasta >> $combined_dir/denv3/combined_denv3.fasta
    elif [[ $serotype == "DENV-4" ]]; then
        cat ${lineage}.fasta >> $combined_dir/denv4/combined_denv4.fasta
    fi
done

echo "Download, renaming, and combination of sequences completed."

# Clean up individual files
rm DENV-*.fasta
