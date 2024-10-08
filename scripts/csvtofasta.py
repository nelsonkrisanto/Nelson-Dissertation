import csv

# Path to your input CSV file
csv_file_path = '/home/people/23203786/scratch/Nelson-Dissertation/primers/cleanned_primers.csv'

# Path to the output FASTA file
fasta_file_path = '/home/people/23203786/scratch/Nelson-Dissertation/primers/dengue_primers.fasta'

# Assuming the primer ID is in the first column, sequence in the second, and genotype in the third
with open(csv_file_path, mode='r') as csv_file, open(fasta_file_path, mode='w') as fasta_file:
    csv_reader = csv.reader(csv_file)
    next(csv_reader, None)  # Skip the header row
    for row in csv_reader:
        primer_id, primer_sequence, genotype = row[0], row[1], row[2]
        fasta_file.write(f'>{primer_id}|{genotype}\n{primer_sequence}\n')

print("Done Converting...")
