#!/bin/bash
#SBATCH --job-name=4.3.3.semi_nested_amp_confusion_matrix
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10

# Load necessary modules
module load python/3.9.15
echo "Python module loaded"

python_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/balance_fasta_sequences_with_metadata.py"
taxonomy_file="/home/people/23203786/scratch/Nelson-Dissertation/results/blast_output/amp_semi_nested/best_taxonomy_assignment_per_sequence.tsv"
fasta_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/unique_passed_semi_nested_sequences.fasta"
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/check/semi"

python "$python_script" "$taxonomy_file" "$fasta_file" "$output_dir"
echo "Python script completed"

# Unload Python module
module unload python/3.9.15
echo "Python module unloaded"

# Load seqkit module
module load seqkit
echo "Seqkit module loaded"

# Split the FASTA file into 10 parts
/home/people/23203786/scratch/Nelson-Dissertation/tools/seqkit split -p 10 "$output_dir/balanced_sequences.fasta"
echo "FASTA file split into 10 parts"

# Unload seqkit module
module unload seqkit
echo "Seqkit module unloaded"

echo "Script completed successfully."

# Load necessary modules
module load python/3.9.15
echo "Python module loaded"
source activate /home/people/23203786/tools/biopython
echo "biopython loaded"

python_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/generate_confusion_matrix.py"
nested_results_folder="/home/people/23203786/scratch/Nelson-Dissertation/results/genome_detective/semi_nested"
genome_detective_fasta_folder="/home/people/23203786/scratch/Nelson-Dissertation/results/genome_detective/semi_nested"
blast_results_file="/home/people/23203786/scratch/Nelson-Dissertation/results/blast_output/amp_semi_nested/best_taxonomy_assignment_per_sequence.tsv"
unique_passed_fasta_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/unique_passed_semi_nested_sequences.fasta"
combined_output_file="/home/people/23203786/scratch/Nelson-Dissertation/results/confusion_matrix/semi_nested/combined_results.csv"
confusion_matrix_output_file="/home/people/23203786/scratch/Nelson-Dissertation/results/confusion_matrix/semi_nested/confusion_matrix.png"

python "$python_script" "$nested_results_folder" "$genome_detective_fasta_folder" "$blast_results_file" "$unique_passed_fasta_file" "$combined_output_file" "$confusion_matrix_output_file"
echo "Python script completed"

source deactivate
module unload python/3.9.15
echo "Python module unloaded"

echo "Script completed successfully."
