# Define the output file path
output_file_path = "/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/protein_db_taxonomy.tsv"

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

taxonomy_template = "Amarillovirales;Flaviviridae;Orthoflavivirus;Orthoflavivirus_denguei;dengue_virus;dengue_virus_type_{};Dengue_virus_{}_{}"

print(f"Writing to output file: {output_file_path}")

with open(output_file_path, "w") as f:
    # Write strain entries
    for strain, accession in dengue_strains.items():
        serotype = ""  # Initialize serotype variable
        if "Dengue virus 1" in strain:
            serotype = "1"
        elif "Dengue virus 2" in strain:
            serotype = "2"
        elif "Dengue virus 3" in strain:
            serotype = "3"
        elif "Dengue virus 4" in strain:
            serotype = "4"
        
        if serotype:
            strain_with_details = strain.replace("Dengue virus ", "").replace(" ", "_")
            taxonomy = taxonomy_template.format(serotype, serotype, strain_with_details)
            formatted_taxonomy = taxonomy.replace(" ", "_")
            line = f"{accession}_DENV{serotype}_{strain_with_details}\t{formatted_taxonomy}\n"
            print(f"Writing line: {line.strip()}")
            f.write(line)
        else:
            print(f"No serotype found for strain: {strain}")

print("Finished writing the output file.")
