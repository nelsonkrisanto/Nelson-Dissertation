#!/bin/bash
#SBATCH --job-name=1.5.create_diamond_db
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nelson.krisanto@ucdconnect.ie
#SBATCH --error=/home/people/23203786/scratch/Nelson-Dissertation/logs/error_%x_%j.txt
#SBATCH --output=/home/people/23203786/scratch/Nelson-Dissertation/logs/log_%x_%j.txt
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# Define the directory for taxonomy files
taxonomy_dir="/home/people/23203786/scratch/Nelson-Dissertation/taxonomy"

# Create the directory if it doesn't exist
mkdir -p $taxonomy_dir

# Create taxon.map
cat <<EOL > $taxonomy_dir/taxon.map
U88536    11053
U87411    11060
M93130    11069
AF326573    11070
EOL

# Create taxon.nodes
cat <<EOL > $taxonomy_dir/taxon.nodes
11053    12637    species
11060    12637    species
11069    12637    species
11070    12637    species
12637    3052464    genus
3052464    1    family
EOL

# Create taxon.names
cat <<EOL > $taxonomy_dir/taxon.names
11053    DENV-1
11060    DENV-2
11069    DENV-3
11070    DENV-4
12637    Dengue virus
3052464    Orthoflavivirus denguei
EOL

echo "Taxonomy files created successfully."
