#!/bin/bash
#SBATCH --job-name=4.2.1.conv_whole_confusion_matrix
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Load necessary modules
module load python/3.9.15
echo "Python module loaded"

# Define paths
amplicons_genotyping_file="/home/people/23203786/scratch/Nelson-Dissertation/results/blast_output/amp/best_taxonomy_assignment_per_sequence.tsv"
whole_genome_genotyping_file="/home/people/23203786/scratch/Nelson-Dissertation/results/blast_output/wholegenome/best_taxonomy_assignment_per_sequence.tsv"
primer_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/filtered_linked_amplicons.csv"
final_csv="/home/people/23203786/scratch/Nelson-Dissertation/results/confusion_matrix/conventional/conv_whole_confusion_matrix.csv"
generate_csv_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/generate_csv_confusion_matrix.py"
confusion_matrix_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/confusion_matrix_whole.py"

# Generate the final CSV for the confusion matrix
echo "Generating the final CSV..."
python $generate_csv_script $amplicons_genotyping_file $whole_genome_genotyping_file $primer_file $final_csv
status=$?
echo "CSV generation script finished with exit code $status"

# Check if the final CSV was generated successfully
if [ $status -ne 0 ] || [ ! -f "$final_csv" ]; then
    echo "Failed to generate the final CSV file: $final_csv"
    exit 1
fi

# Print the contents of the final CSV for debugging
echo "Contents of the final CSV file:"
cat $final_csv

# Run the confusion matrix script
echo "Running confusion matrix script..."
python $confusion_matrix_script $final_csv
status=$?
echo "Confusion matrix script finished with exit code $status"

# Unload modules
module unload python/3.9.15
echo "Python module unloaded"

echo "Confusion matrix generation completed successfully."
