#!/bin/bash
#SBATCH --job-name=1.4.3.extract_vadr_amplicons_semi_nested
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Define the output directory and the pattern to search for .sqa and .sqc files
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/vadr_genotyping_insilico_semi_nested/vadr_output"
sqa_pattern="$output_dir/*.vadr.sqa"
sqc_pattern="$output_dir/*.vadr.sqc"

# Define the combined output files
combined_sqa_file="$output_dir/vadr_combined.sqa"
combined_sqc_file="$output_dir/vadr_combined.sqc"

# Merge all .sqa files
cat $sqa_pattern > $combined_sqa_file

# Merge all .sqc files
cat $sqc_pattern > $combined_sqc_file

# Verify the merge
if [[ -f "$combined_sqa_file" && -f "$combined_sqc_file" ]]; then
    echo "Merged .sqa and .sqc files created successfully."
else
    echo "Failed to create the merged .sqa and .sqc files."
    exit 1
fi

# Load Python module
module load python/3.9.15

# Set the paths for the input and output files
vadr_sqa_file="/home/people/23203786/scratch/Nelson-Dissertation/results/vadr_genotyping_insilico_semi_nested/vadr_output/vadr_combined.sqa"
vadr_sqc_file="/home/people/23203786/scratch/Nelson-Dissertation/results/vadr_genotyping_insilico_semi_nested/vadr_output/vadr_combined.sqc"
original_fasta="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/filtered_semi_nested_insilico_PCR_amplicons.fasta"
passed_fasta="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/passed_semi_nested_sequences.fasta"
python_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/extract_vadr.py"

# Run the Python script to extract passed sequences
python $python_script $vadr_sqa_file $vadr_sqc_file $original_fasta $passed_fasta

# Unload Python module
module unload python/3.9.15