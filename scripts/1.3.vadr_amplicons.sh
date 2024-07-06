#!/bin/bash
#SBATCH --job-name=1.3.vadr_amplicons
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

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
