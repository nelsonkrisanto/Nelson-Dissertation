#!/bin/bash
#SBATCH --job-name=02.fetch_papers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=5

# Directories
my_dir="/home/people/23203786/scratch/Nelson-Dissertation"
reference_path="/home/people/23203786/scratch/Nelson-Dissertation/references"
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/primers"

# Load necessary modules
module load bwa/0.7
module load samtools/1.10

# Change to the output directory
cd "$output_dir"
if [ $? -ne 0 ]; then
    echo "Failed to change directory to $output_dir"
    exit 1
fi

# Define file paths
primers_file="$my_dir/primers/dengue_primers.fasta"
output_tsv="mapping_positions.tsv"

# Create the header for the TSV file
echo -e "Primer\tReference\tStart\tEnd" > "$output_tsv"
if [ $? -ne 0 ]; then
    echo "Failed to write header to $output_tsv"
    exit 1
fi

# Loop to read each primer and reference and perform the alignment
while IFS= read -r line
do
    if [[ "$line" =~ ^\> ]]; then
        primer_id=$(echo "$line" | tr -d '>' | tr -d '\r')  # Remove '>' and carriage return
    else
        primer_sequence="$line"

        # Align each primer to each reference genome
        for ref_file in "$reference_path"/*.fna; do
            ref_name=$(basename "$ref_file" .fna)
            
            # Create temporary FASTA file for primer
            primer_fasta="$output_dir/temp_primer_${primer_id}.fasta"
            echo -e ">$primer_id\n$primer_sequence" > "$primer_fasta"

            # Perform BWA alignment and log any errors
            sai_file="$output_dir/${primer_id}_${ref_name}.sai"
            bwa aln -o 0 -n 9 -l 100 -N "$ref_file" "$primer_fasta" > "$sai_file" 2>> "$output_dir/bwa_error_${primer_id}_${ref_name}.log"
            if [ $? -ne 0 ]; then
                echo "BWA alignment failed for primer $primer_id and reference $ref_name"
                exit 1
            fi

            # ... (rest of the alignment and conversion to SAM)
        done
    fi
done < "$primers_file"

echo "Mapping positions extracted and saved to $output_tsv."


# Unload modules
module unload bwa/0.7
module unload samtools/1.10
