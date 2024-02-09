#!/bin/bash

# Set project directories and paths
my_dir="/home/people/23203786/scratch/_temp"
capsid_path="$my_dir/raw_data/reference_sequences_VP1"
rdrp_path="$my_dir/raw_data/reference_sequences_RdRp"
scripts_path="$my_dir/scripts"
scrappaper="/home/people/23203786/tools/ScrapPaper/scrappaper.py"
entrezdirect="/home/people/23203786/tools/entrez-direct_env"

# Change working directory to the project directory
cd "$my_dir"

# Ensure project directory structure
project_structure=("data" "logs" "raw_data" "results")

# Create subdirectories if they don't exist
for subdir in "${project_structure[@]}"; do
    if [ ! -d "$subdir" ]; then
        mkdir -p "$subdir"
    fi
done

# Create subdirectories for reference sequences
if [ ! -d "$capsid_path" ]; then
    mkdir -p "$capsid_path"
fi

if [ ! -d "$rdrp_path" ]; then
    mkdir -p "$rdrp_path"
fi

# Print message indicating successful setup
echo "Project subdirectories are set up in $my_dir."

# Check if the scripts directory exists
if [ -d "$scripts_path" ]; then
    # Print message indicating that the scripts directory already exists
    echo "Scripts directory $scripts_path already exists."
else
    # Create scripts directory if it doesn't exist
    mkdir -p "$scripts_path"

    # Print message indicating successful setup
    echo "Scripts directory is set up in $scripts_path."
fi