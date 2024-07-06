#!/bin/bash
#SBATCH --job-name=1.4.1.extract_vadr_cls
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Define the base output directory
base_output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/vadr_genotyping_cls/vadr_output"

# Define the combined output files
combined_sqa_file="$base_output_dir/vadr_combined.sqa"
combined_sqc_file="$base_output_dir/vadr_combined.sqc"

# Initialize combined files
> "$combined_sqa_file"
> "$combined_sqc_file"

# Loop through each vadr_output folder and concatenate .sqa and .sqc files
for dir in "$base_output_dir"/*; do
    if [ -d "$dir" ]; then
        sqa_file="$dir/$(basename "$dir").vadr.sqa"
        sqc_file="$dir/$(basename "$dir").vadr.sqc"

        if [ -f "$sqa_file" ]; then
            cat "$sqa_file" >> "$combined_sqa_file"
        fi

        if [ -f "$sqc_file" ]; then
            cat "$sqc_file" >> "$combined_sqc_file"
        fi
    fi
done

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
vadr_sqa_file="$base_output_dir/vadr_combined.sqa"
vadr_sqc_file="$base_output_dir/vadr_combined.sqc"
original_fasta="/home/people/23203786/scratch/Nelson-Dissertation/results/mmseqs2_clustering/clustered_sequences.fasta"
passed_fasta="/home/people/23203786/scratch/Nelson-Dissertation/results/wholegenome/passed_whole_genome_sequences.fasta"
python_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/extract_vadr.py"

# Run the Python script to extract passed sequences
python $python_script $vadr_sqa_file $vadr_sqc_file $original_fasta $passed_fasta

# Unload Python module
module unload python/3.9.15
