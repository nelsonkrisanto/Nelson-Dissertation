# Revised script to process fasta files and update headers

import os

# Dictionary with strain names and their corresponding protein accession numbers
dengue_strains = {
    "Dengue virus 1 Brazil/97-11/1997": "P27909",
    "Dengue virus 1 Jamaica/CV1636/1977": "P27913",
    "Dengue virus 1 Nauru/West Pac/1974": "P17763",
    "Dengue virus 1 Singapore/S275/1990": "P33478",
    "Dengue virus 1 Thailand/AHF 82-80/1980": "P27912",
    "Dengue virus 2 16681-PDK53": "P29991",
    "Dengue virus 2 China/D2-04": "P30026",
    "Dengue virus 2 Jamaica/1409/1983": "P07564",
    "Dengue virus 2 Malaysia M2": "P14338",
    "Dengue virus 2 Malaysia M3": "P14339",
    "Dengue virus 2 Peru/IQT2913/1996": "Q9WDA6",
    "Dengue virus 2 Puerto Rico/PR159-S1/1969": "P12823",
    "Dengue virus 2 Thailand/0168/1979": "P14337",
    "Dengue virus 2 Thailand/16681/84": "P29990",
    "Dengue virus 2 Thailand/NGS-C/1944": "P14340",
    "Dengue virus 2 Thailand/PUO-218/1980": "P18356",
    "Dengue virus 2 Thailand/TH-36/1958": "P29984",
    "Dengue virus 2 Tonga/EKB194/1974": "P27914",
    "Dengue virus 3 China/80-2/1980": "Q99D35",
    "Dengue virus 3 Martinique/1243/1999": "Q6YMS3",
    "Dengue virus 3 Philippines/H87/1956": "P27915",
    "Dengue virus 3 Singapore/8120/1995": "Q5UB51",
    "Dengue virus 3 Sri Lanka/1266/2000": "Q6YMS4",
    "Dengue virus 4 Dominica/814669/1981": "P09866",
    "Dengue virus 4 Philippines/H241/1956": "Q58HT7",
    "Dengue virus 4 Singapore/8976/1995": "Q5UCB8",
    "Dengue virus 4 Thailand/0348/1991": "Q2YHF0",
    "Dengue virus 4 Thailand/0476/1997": "Q2YHF2"
}

taxonomy_template = "Amarillovirales;Flaviviridae;Orthoflavivirus;Orthoflavivirus denguei;dengue virus;dengue virus type {};{}"

# Define output directories for each serotype
output_dirs = {
    "denv1": "/home/people/23203786/scratch/Nelson-Dissertation/db/denv1/combined_denv1.fasta",
    "denv2": "/home/people/23203786/scratch/Nelson-Dissertation/db/denv2/combined_denv2.fasta",
    "denv3": "/home/people/23203786/scratch/Nelson-Dissertation/db/denv3/combined_denv3.fasta",
    "denv4": "/home/people/23203786/scratch/Nelson-Dissertation/db/denv4/combined_denv4.fasta"
}

# Combine fasta files and update headers
for serotype, output_file in output_dirs.items():
    with open(output_file, "w") as outfile:
        for strain, accession in dengue_strains.items():
            if serotype in accession.lower():
                fasta_file = f"/home/people/23203786/scratch/Nelson-Dissertation/db/{serotype}/{accession}.fasta"
                if os.path.exists(fasta_file):
                    with open(fasta_file, "r") as infile:
                        for line in infile:
                            if line.startswith(">"):
                                strain_name = strain.replace(" ", "_").replace("/", "_")
                                header = f">{accession}_{serotype.upper()}_{strain_name}\n"
                                outfile.write(header)
                            else:
                                outfile.write(line)
                else:
                    print(f"Warning: Fasta file {fasta_file} not found.")

print("Fasta files combined and headers updated.")
