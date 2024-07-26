#!/bin/bash
#SBATCH --job-name=allconv
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Load the Perl module
module load perl/5.30 

# Define input file paths
primers="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/conventional_primer_combinations.tsv"
sequences="/home/people/23203786/scratch/Nelson-Dissertation/results/fasta/combined_sequences.fasta"
out_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico"
tool="/home/people/23203786/tools/in_silico_PCR/in_silico_PCR.pl"
filter_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/filter_amplicons.py"

# Create output directory if it doesn't exist
mkdir -p "$out_dir"

# Change working directory
cd "$out_dir" || exit

# Make sure the Perl script is executable
chmod +x "$tool"

# Run in-silico PCR with options to exclude primer sequences, allow 1-2 mismatches for degenerate bases, and ensure correct orientation
perl "$tool" -s "$sequences" -p "$primers" -l 1200 -m -i -e -r > conv_insilico_PCR_results.txt 2> conv_insilico_PCR_amplicons.fasta

# Unload Perl module
module unload perl/5.30

# Load Python module
module load python/3.9.15

# Run the filtering script to filter amplicons between 300-500bp
python "$filter_script" conv_insilico_PCR_results.txt conv_insilico_PCR_amplicons.fasta filtered_conv_insilico_PCR_amplicons.fasta

# Unload Python module
module unload python/3.9.15

# Load the Anaconda module
module load anaconda/3.5.2

# Activate the Conda environment for VADR
source /opt/software/anaconda/3.5.2/etc/profile.d/conda.sh
conda activate vadr

# Set PERL5LIB to include the directories of vadr.pm and sqp_opts.pm explicitly
export PERL5LIB=$CONDA_PREFIX/share/vadr-1.5.1/vadr:$CONDA_PREFIX/share/sequip-0.10/lib:$PERL5LIB

# Verify PERL5LIB paths
echo "PERL5LIB is set to: $PERL5LIB"

# Check if the v-annotate.pl script is available
vadr_script=$(find $CONDA_PREFIX -name v-annotate.pl | head -n 1)

if [[ -z "$vadr_script" ]]; then
    echo "v-annotate.pl script not found. Please check the VADR installation."
    conda deactivate
    exit 1
fi

# Check if vadr.pm is available
vadr_pm=$(find $CONDA_PREFIX -name vadr.pm | head -n 1)

if [[ -z "$vadr_pm" ]]; then
    echo "vadr.pm not found. Please check the VADR installation."
    conda deactivate
    exit 1
else
    echo "Found vadr.pm at $vadr_pm"
fi

# Check if sqp_opts.pm is available
sqp_opts_pm=$(find $CONDA_PREFIX -name sqp_opts.pm | head -n 1)

if [[ -z "$sqp_opts_pm" ]]; then
    echo "sqp_opts.pm not found. Please check the SEQUIP installation."
    conda deactivate
    exit 1
else
    echo "Found sqp_opts.pm at $sqp_opts_pm"
fi

# Define input and output file paths
input_fasta="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/filtered_conv_insilico_PCR_amplicons.fasta"
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/vadr_genotyping_insilico"
mkdir -p "$output_dir"

# Define the correct VADR model directory and output prefix
vadr_model_dir="/home/people/23203786/.conda/envs/vadr/share/vadr/vadr-models/vadr-models-flavi-1.2-1"
output_prefix="$output_dir/vadr_output"

# Check if the model directory exists and contains the required files
if [ ! -d "$vadr_model_dir" ] || [ ! -f "$vadr_model_dir/flavi.minfo" ]; then
    echo "Model directory or required files not found in $vadr_model_dir. Please check the VADR installation."
    conda deactivate
    exit 1
else
    echo "Model directory and required files found in $vadr_model_dir"
fi

# Run VADR with Flaviviridae models for Dengue virus
perl -I$CONDA_PREFIX/share/vadr-1.5.1/vadr -I$CONDA_PREFIX/share/sequip-0.10/lib $vadr_script --split --cpu 10 --group Dengue --nomisc --noprotid --mkey flavi --mdir $vadr_model_dir --f $input_fasta $output_prefix

# Deactivate Conda environment
conda deactivate

# Define the output directory and the pattern to search for .sqa and .sqc files
output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/vadr_genotyping_insilico/vadr_output"
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
vadr_sqa_file="/home/people/23203786/scratch/Nelson-Dissertation/results/vadr_genotyping_insilico/vadr_output/vadr_combined.sqa"
vadr_sqc_file="/home/people/23203786/scratch/Nelson-Dissertation/results/vadr_genotyping_insilico/vadr_output/vadr_combined.sqc"
original_fasta="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/filtered_conv_insilico_PCR_amplicons.fasta"
passed_fasta="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/passed_sequences.fasta"
python_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/extract_vadr.py"

# Run the Python script to extract passed sequences
python $python_script $vadr_sqa_file $vadr_sqc_file $original_fasta $passed_fasta

# Unload Python module
module unload python/3.9.15

# Load necessary modules
module load blast/2.12.0
source activate /home/people/23203786/tools/biopython

# Define paths
query_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/passed_sequences.fasta"
unique_query_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/unique_passed_sequences.fasta"  # Path for unique sequences
base_output_dir="/home/people/23203786/scratch/Nelson-Dissertation/results/blast_output/amp"
taxonomy_script="/home/people/23203786/scratch/Assign-Taxonomy-with-BLAST/taxonomy_assignment_BLAST.py"
taxonomy_file="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/protein_db_taxonomy.tsv"
best_taxonomy_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/best_taxonomy_assignment.py"
seqkit_path="/home/people/23203786/scratch/Nelson-Dissertation/tools/seqkit"

# Create a temporary directory and set it as TMPDIR
TMP_DIR=$(mktemp -d)
export TMPDIR=$TMP_DIR
trap 'rm -rf "$TMP_DIR"' EXIT

# Ensure the base output directory exists
mkdir -p $base_output_dir

# Remove duplicate sequences
echo "Removing duplicate sequences..."
$seqkit_path rmdup -s -i $query_file -o $unique_query_file

# Define serotypes and their database paths
declare -A serotypes=(
    ["denv4"]="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/denv4/dengue_blast_db"
    ["denv3"]="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/denv3/dengue_blast_db"
    ["denv2"]="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/denv2/dengue_blast_db"
    ["denv1"]="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/RDP_output/denv1/dengue_blast_db"
)

# Run BLAST and taxonomy assignment for each serotype
for serotype in "${!serotypes[@]}"; do
    db_path="${serotypes[$serotype]}"
    blast_output="${base_output_dir}/${serotype}_blast_output.txt"
    
    # Run BLAST search
    echo "Running BLAST search for $serotype..."
    blastx -db $db_path -query $unique_query_file -out $blast_output -evalue 1e-5 -outfmt "6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore staxids" -num_threads 10

    # Check if BLAST output file is generated
    if [ ! -f "$blast_output" ]; then
        echo "BLAST search for $serotype failed. No output file generated."
        exit 1
    else
        echo "BLAST search for $serotype completed successfully. Output file: $blast_output"
    fi

    # Run taxonomy assignment
    echo "Running taxonomy assignment for $serotype..."
    python3 $taxonomy_script --blast_flavor blastx --cutoff_species 70 --cutoff_family 60 --cutoff_phylum 50 --output_dir $base_output_dir --verbose --blast_file $blast_output -- $unique_query_file $taxonomy_file
    if [ $? -ne 0 ]; then
        echo "Taxonomy assignment for $serotype failed."
        exit 1
    else
        echo "Taxonomy assignment for $serotype completed successfully."
    fi
done

# Combine taxonomy assignments and determine the best one for each sequence
echo "Determining the best taxonomy assignment for each sequence..."
combined_taxonomy_file="${base_output_dir}/taxonomy_assignment_per_sequence.tsv"
best_output_file="${base_output_dir}/best_taxonomy_assignment_per_sequence.tsv"

python3 $best_taxonomy_script $combined_taxonomy_file $best_output_file

# Unload modules
module unload blast/2.12.0
source deactivate

echo "All BLAST searches and taxonomy assignments completed successfully."
echo "Best taxonomy assignment per sequence has been saved to ${best_output_file}"

# Load Python module
module load python/3.9.15

# Define input file paths
results_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/conv_insilico_PCR_results.txt"
amplicons_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/conv_insilico_PCR_amplicons.fasta"
primers_file="/home/people/23203786/scratch/Nelson-Dissertation/results/tsv/conventional_primer_combinations.tsv"
output_file="/home/people/23203786/scratch/Nelson-Dissertation/results/insilico/filtered_linked_amplicons.csv"

# Path to the link_amplicons script
link_script="/home/people/23203786/scratch/Nelson-Dissertation/scripts/link_conv_amplicons.py"

# Run the linking script to filter amplicons between 300-1200bp and link them to their accession numbers and primer pairs
python "$link_script" "$results_file" "$amplicons_file" "$primers_file" "$output_file"

# Unload Python module
module unload python/3.9.15
