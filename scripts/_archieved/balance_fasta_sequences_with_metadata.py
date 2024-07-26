import os
import random
from collections import defaultdict

def read_taxonomy_file(taxonomy_file):
    sequences_by_genotype = defaultdict(list)
    taxonomy_data = {}
    with open(taxonomy_file, 'r') as file:
        next(file)  # Skip header line if it exists
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                amp_id = parts[0]
                taxonomy = parts[1].split(';')
                genotype = taxonomy[5]  # dengue_virus_type_X
                serotype = taxonomy[6] if len(taxonomy) > 6 else 'undetermined'
                sequences_by_genotype[genotype].append(amp_id)
                taxonomy_data[amp_id] = (genotype, serotype)
    return sequences_by_genotype, taxonomy_data

def select_sequences(sequences_by_genotype, max_sequences_per_genotype=500):
    selected_sequences = []
    for genotype, sequences in sequences_by_genotype.items():
        if len(sequences) > max_sequences_per_genotype:
            selected_sequences.extend(random.sample(sequences, max_sequences_per_genotype))
        else:
            selected_sequences.extend(sequences)
    return selected_sequences

def read_fasta(fasta_file):
    sequences = {}
    with open(fasta_file, 'r') as file:
        header = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequences[header] = ''.join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if header:
            sequences[header] = ''.join(sequence)
    return sequences

def write_fasta(sequences, selected_ids, output_file):
    with open(output_file, 'w') as file:
        for seq_id in selected_ids:
            if seq_id in sequences:
                file.write(f">{seq_id}\n")
                file.write(f"{sequences[seq_id]}\n")

def write_metadata(selected_ids, taxonomy_data, metadata_file):
    with open(metadata_file, 'w') as file:
        file.write("Amplicons\tGenotype\tSerotype\n")
        for seq_id in selected_ids:
            if seq_id in taxonomy_data:
                genotype, serotype = taxonomy_data[seq_id]
                file.write(f"{seq_id}\t{genotype}\t{serotype}\n")

def main(taxonomy_file, fasta_file, output_dir, max_sequences_per_genotype=500):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sequences_by_genotype, taxonomy_data = read_taxonomy_file(taxonomy_file)
    selected_sequences = select_sequences(sequences_by_genotype, max_sequences_per_genotype)
    sequences = read_fasta(fasta_file)
    
    output_fasta_file = os.path.join(output_dir, 'balanced_sequences.fasta')
    write_fasta(sequences, selected_sequences, output_fasta_file)
    
    output_metadata_file = os.path.join(output_dir, 'balanced_sequences_metadata.tsv')
    write_metadata(selected_sequences, taxonomy_data, output_metadata_file)

    print(f"Balanced sequences written to {output_fasta_file}")
    print(f"Metadata written to {output_metadata_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python balance_fasta_sequences_with_metadata.py taxonomy_file fasta_file output_dir")
        sys.exit(1)
    taxonomy_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_dir = sys.argv[3]
    main(taxonomy_file, fasta_file, output_dir)
