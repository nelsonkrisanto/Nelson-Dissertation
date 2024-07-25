#!/bin/bash
#SBATCH --job-name=4.3.1.conv_amp_confusion_matrix
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10

# Load necessary modules
module load python/3.9.15
echo "Python module loaded"
source activate /home/people/23203786/tools/biopython
echo "biopython loaded"

python_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/generate_confusion_matrix.py"
nested_results_folder="/home/people/23203786/scratch/Nelson-Dissertation/results/genome_detective/conv"
genome_detective_fasta_folder="/home/people/23203786/scratch/Nelson-Dissertation/results/genome_detective/conv"
blast_results_file="/home/people/23203786/scratch/Nelson-Dissertation/results/blast_output/amp/best_taxonomy_assignment_per_sequence.tsv"
unique_passed_fasta_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/unique_passed_sequences.fasta"
combined_output_file="/home/people/23203786/scratch/Nelson-Dissertation/results/confusion_matrix/conventional/combined_results.csv"
confusion_matrix_output_file="/home/people/23203786/scratch/Nelson-Dissertation/results/confusion_matrix/conventional/confusion_matrix.png"

python "$python_script" "$nested_results_folder" "$genome_detective_fasta_folder" "$blast_results_file" "$unique_passed_fasta_file" "$combined_output_file" "$confusion_matrix_output_file"
echo "Python script completed"

source deactivate
module unload python/3.9.15
echo "Python module unloaded"

echo "Script completed successfully."