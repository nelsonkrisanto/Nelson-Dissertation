#!/bin/bash
#SBATCH --job-name=blasttax
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Define paths
taxonomy_dir="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy"
mkdir -p $taxonomy_dir

# Generate .tax files for each serotype
for serotype in denv1 denv2 denv3 denv4; do
    fasta_file="/home/people/23203786/scratch/Nelson-Dissertation/db/$serotype/$serotype.fasta"
    tax_file="$taxonomy_dir/$serotype.tax"
    
    # Determine the taxonomic ID
    case $serotype in
        denv1)
            taxid=11053 ;;
        denv2)
            taxid=11060 ;;
        denv3)
            taxid=11069 ;;
        denv4)
            taxid=11070 ;;
    esac
    
    # Generate .tax file
    awk -v taxid=$taxid '/^>/ {print substr($1, 2) " " taxid}' $fasta_file > $tax_file
done

echo "Taxonomy files created successfully."
