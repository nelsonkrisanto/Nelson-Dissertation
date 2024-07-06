#!/bin/bash
#SBATCH --job-name=1.8.1.whole_genome_confusion_matrix
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
file="/home/people/23203786/scratch/Nelson-Dissertation/results/confusion_matrix/whole_genome_confusion.csv"
python_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/whole_genome_confusion_matrix.py"
echo "Paths defined: file=$file, python_script=$python_script"

# Check if the file exists
if [ ! -f "$file" ]; then
    echo "CSV file not found: $file"
    exit 1
fi

# Run the Python script
echo "Running Python script..."
python $python_script $file
status=$?
echo "Python script finished with exit code $status"

# Unload modules
module unload python/3.9.15
echo "Python module unloaded"

echo "Confusion matrix generation completed successfully."
