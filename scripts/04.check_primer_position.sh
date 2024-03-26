#!/bin/bash
#SBATCH --job-name=04.check_primer_position
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

# Ensure the output directory exists
mkdir -p "$output_dir" || { echo "Failed to create output directory"; exit 1; }

# Change to the output directory
cd "$output_dir" || { echo "Failed to change directory to $output_dir"; exit 1; }

# Define file paths
primers_file="$my_dir/primers/dengue_primers.fasta"
output_tsv="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/mapping_positions.tsv"

# Create the header for the TSV file
echo -e "Primer\tReference\tStart\tEnd" > "$output_tsv" || { echo "Failed to write header to $output_tsv"; exit 1; }

# Loop to read each primer and reference and perform the alignment
while IFS= read -r line
do
    if [[ "$line" =~ ^\> ]]; then
        # Remove '>' and potential carriage returns, and replace '/' with '_'
        primer_id=$(echo "$line" | tr -d '>' | tr -d '\r' | tr '/' '_')
        primer_sequence=""
    else
        primer_sequence="$line"
        
        # Align each primer to each reference genome
        for ref_file in "$reference_path"/*.fna; do
            ref_name=$(basename "$ref_file" .fna)
            
            # Check if index exists, if not, create it
            if [ ! -f "${ref_file}.bwt" ]; then
                bwa index "$ref_file" || { echo "Indexing failed for $ref_file"; exit 1; }
            fi
            
            # Create temporary FASTA file for primer
            primer_fasta="$output_dir/temp_primer_${primer_id}.fasta"
            echo -e ">$primer_id\n$primer_sequence" > "$primer_fasta" || { echo "Failed to create FASTA for $primer_id"; exit 1; }

            # Perform BWA alignment and log any errors
            sai_file="$output_dir/${primer_id}_${ref_name}.sai"
            bwa aln -o 0 -n 9 -l 100 -N "$ref_file" "$primer_fasta" > "$sai_file" 2>> "$output_dir/bwa_error_${primer_id}_${ref_name}.log"
            if [ $? -ne 0 ]; then
                echo "BWA alignment failed for primer $primer_id and reference $ref_name"
                exit 1
            fi

            # Convert to SAM format and extract the alignment info
            sam_file="$output_dir/${primer_id}_${ref_name}.sam"
            bwa samse "$ref_file" "$sai_file" "$primer_fasta" > "$sam_file" 2>> "$output_dir/bwa_error_${primer_id}_${ref_name}.log"
            if [ $? -ne 0 ]; then
                echo "BWA samse failed for primer $primer_id and reference $ref_name"
                exit 1
            fi

            # Extract alignment positions and append to the TSV file, omitting if start position is 0
            alignment_info=$(grep -v "^@" "$sam_file" | awk -v ref="$ref_name" -v primer="$primer_id" '$4 > 0 {print primer "\t" ref "\t" $4 "\t" $4 + length($10) - 1}')
            
            # Check if alignment_info is not empty before appending to the TSV
            if [ -n "$alignment_info" ]; then
                echo -e "$alignment_info" >> "$output_tsv"
            fi

            # Clean up the temporary SAM file
            rm "$sam_file"
        done
    fi
done < "$primers_file"

echo "Mapping positions extracted and saved to $output_tsv."

# Unload modules
module unload bwa/0.7
module unload samtools/1.10
