# Define the output file path
output_file_path = "/home/people/23203786/scratch/Nelson-Dissertation/taxonomy/protein_db_taxonomy.tsv"

# Dictionary with strain names and their corresponding protein accession numbers
dengue_strains = {
    "DENV-1_genotype_V": "P27909",
    "DENV-1_genotype_II": "P27913",
    "DENV-1_genotype_I": "P17763",
    "DENV-1_genotype_IV": "P33478",
    "DENV-1_genotype_III": "P27912",
    "DENV-2_genotype_Asian_I": "P29991",
    "DENV-2_genotype_Asian_II": "P30026",
    "DENV-2_genotype_American": "P07564",
    "DENV-2_genotype_Cosmopolitan": "P14338",
    "DENV-3_genotype_III": "Q99D35",
    "DENV-3_genotype_IV": "Q6YMS3",
    "DENV-3_genotype_I": "P27915",
    "DENV-3_genotype_II": "Q5UB51",
    "DENV-3_genotype_V": "Q6YMS4",
    "DENV-4_genotype_II": "P09866",
    "DENV-4_genotype_I": "Q58HT7"
}

taxonomy_template = "Amarillovirales;Flaviviridae;Orthoflavivirus;Orthoflavivirus_denguei;dengue_virus;dengue_virus_type_{};{}"

print(f"Writing to output file: {output_file_path}")

with open(output_file_path, "w") as f:
    # Write strain entries
    for strain, accession in dengue_strains.items():
        serotype = strain.split("_")[0].split("-")[1]  # Extract the serotype from the strain name
        taxonomy = taxonomy_template.format(serotype, strain)
        line = f"{strain}\t{taxonomy}\n"
        print(f"Writing line: {line.strip()}")
        f.write(line)

print("Finished writing the output file.")
