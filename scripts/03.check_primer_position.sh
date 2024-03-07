#!/bin/bash
#SBATCH --job-name=02.fetch_papers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie,nelson_krisanto@hotmail.com # Where to send mail
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt  # Error log
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt  # Standard output log
#SBATCH --cpus-per-task=5

module load bwa/0.7
module load samtools/1.10

cd "$output_dir"

primers_file="primers.fasta"   # Replace with your primer sequences FASTA file
output_tsv="mapping_positions.tsv"  # Output TSV file name

# Create output directory
mkdir -p "$output_dir"

# Create the header for the TSV file
echo -e "Primer\tReference\tStart\tEnd" > "$output_tsv"

# Read each line in the primers file
while IFS= read -r line
do
    if [[ "$line" =~ ^\> ]]; then
        primer_id=$(echo "$line" | tr -d '>')
        primer_sequence=""
    else
        primer_sequence="$line"

        for ref_file in "$REFERENCE_PATH"/*.fasta; do
            ref_name=$(basename "$ref_file" .fasta)
            
            # Index the reference genome if not already indexed
            bwa index "$ref_file"

            # Perform BWA alignment with more sensitive settings
            primer_fasta="$output_dir/temp_primer_${primer_id}.fasta"
            echo -e ">$primer_id\n$primer_sequence" > "$primer_fasta"
            bwa aln -o 0 -n 9 -l 100 -N "$ref_file" "$primer_fasta" > "$output_dir/${primer_id}_${ref_name}.sai"
            bwa samse "$ref_file" "$output_dir/${primer_id}_${ref_name}.sai" "$primer_fasta" > "$output_dir/${primer_id}_${ref_name}.sam"
            rm "$primer_fasta" "$output_dir/${primer_id}_${ref_name}.sai"

            # Extract alignment positions and append to the TSV file, omitting if start position is 0
            alignment_info=$(grep -v "^@" "$output_dir/${primer_id}_${ref_name}.sam" | awk -v ref="$ref_name" '$4 > 1000 {print $4 "\t" $4 + length($10) - 1}')
            
            # Check if alignment_info is not empty before appending to the TSV
            if [ -n "$alignment_info" ]; then
                echo -e "$primer_id\t$ref_name\t$alignment_info" >> "$output_tsv"
            fi
        done
    fi
done < "$primers_file"

echo "Mapping positions extracted and saved to $output_tsv."

module unload bwa/0.7
module unload samtools/1.10
