#!/bin/bash
#SBATCH --job-name=convert_to_RDP
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00

# Load Python module
module load python/3.7.4

# Set PYTHONPATH to include FROGS library
export PYTHONPATH=/home/people/23203786/scratch/FROGS/lib:$PYTHONPATH

# Define paths
python_script="/home/people/23203786/scratch/FROGS/libexec/fasta2RDP.py"
input_fasta="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/refgenome.fasta"
input_taxonomy="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/protein_db_taxonomy.tsv"
output_rdp_fasta="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/protein_db_taxonomy_RDP.fasta"
output_rdp_taxonomy="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/protein_db_taxonomy_RDP.tax"

# Run the Python script
python $python_script -d $input_fasta -t $input_taxonomy -r Domain Phylum Class Order Family Genus Species --rdp-taxonomy $output_rdp_taxonomy --rdp-fasta $output_rdp_fasta

echo "Conversion to RDP format completed."
