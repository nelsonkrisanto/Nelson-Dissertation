#!/bin/bash
#SBATCH --job-name=3.2.csvtofasta
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=2

# Load Python module
module load python/3.9.15

# Run the Python script
python /home/people/23203786/scratch/Nelson-Dissertation/scripts/3.2.csvtofasta.py

# Unload Python module
module unload python/3.9.15
