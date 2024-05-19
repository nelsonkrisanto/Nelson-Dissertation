#!/bin/bash
#SBATCH --job-name=3.4.fetch_temp
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=2

# Load Anaconda module
module load anaconda/3.5.2

# Activate the conda environment for Biopython
source activate /home/people/23203786/tools/biopython

# Ensure necessary packages are installed
conda install -y pandas biopython tqdm openpyxl

# Run the Python script
python /home/people/23203786/scratch/Nelson-Dissertation/scripts/fetch_temp.py

# Deactivate the conda environment
conda deactivate

# Unload Anaconda module
module unload anaconda/3.5.2
