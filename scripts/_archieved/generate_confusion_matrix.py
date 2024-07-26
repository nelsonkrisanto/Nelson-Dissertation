import matplotlib
matplotlib.use('Agg')

import pandas as pd
import glob
import os
import sys
from Bio import SeqIO
from sklearn.metrics import confusion_matrix, classification_report
import matplotlib.pyplot as plt
import seaborn as sns

# Complete mapping from Genome Detective to new lineage nomenclature
genotype_mapping = {
    '2V': 'DENV-2_genotype_V',
    '2III': 'DENV-2_genotype_III',
    '2II': 'DENV-2_genotype_Asian_II',
    '2I': 'DENV-2_genotype_Asian_I',
    '2AM': 'DENV-2_genotype_American',
    '2CO': 'DENV-2_genotype_Cosmopolitan',
    '1V': 'DENV-1_genotype_V',
    '1II': 'DENV-1_genotype_II',
    '1I': 'DENV-1_genotype_I',
    '1IV': 'DENV-1_genotype_IV',
    '1III': 'DENV-1_genotype_III',
    '3V': 'DENV-3_genotype_V',
    '3IV': 'DENV-3_genotype_IV',
    '3III': 'DENV-3_genotype_III',
    '3II': 'DENV-3_genotype_II',
    '3I': 'DENV-3_genotype_I',
    '4II': 'DENV-4_genotype_II',
    '4I': 'DENV-4_genotype_I'
}

def extract_serotype(species_column):
    species_column = species_column.lower().replace('_', ' ').replace('-', ' ')
    if 'dengue virus type 1' in species_column or 'dengue virus type i' in species_column:
        return 'DENV1'
    elif 'dengue virus type 2' in species_column or 'dengue virus type ii' in species_column:
        return 'DENV2'
    elif 'dengue virus type 3' in species_column or 'dengue virus type iii' in species_column:
        return 'DENV3'
    elif 'dengue virus type 4' in species_column or 'dengue virus type iv' in species_column:
        return 'DENV4'
    return 'Unknown'

def map_genotype(assignment_column):
    genotype_key = assignment_column.split('_')[0]
    mapped_genotype = genotype_mapping.get(genotype_key, 'Unknown')
    print(f"Mapping genotype: {genotype_key} to {mapped_genotype}")
    return mapped_genotype

def read_fasta_sequences(fasta_file):
    fasta_sequences = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        seq = str(record.seq).replace("\n", "").strip().lower()
        fasta_sequences[record.id.strip()] = seq
    return fasta_sequences

def generate_confusion_matrix(nested_results_folder, genome_detective_fasta_folder, blast_results_file, unique_passed_fasta_file, combined_output_file, confusion_matrix_output_file):
    print(f"Nested results folder: {nested_results_folder}")
    print(f"Genome detective FASTA folder: {genome_detective_fasta_folder}")
    print(f"BLAST results file: {blast_results_file}")
    print(f"Unique passed FASTA file: {unique_passed_fasta_file}")
    print(f"Combined output file: {combined_output_file}")
    print(f"Confusion matrix output file: {confusion_matrix_output_file}")

    # Step 1: Read all CSV files from the nested results folder
    nested_files = glob.glob(os.path.join(nested_results_folder, '*.csv'))
    nested_dfs = []
    for file in nested_files:
        try:
            df = pd.read_csv(file, engine='python', on_bad_lines='warn')
            nested_dfs.append(df)
            print(f"Loaded {file} with {len(df)} entries")
        except Exception as e:
            print(f"Failed to read {file}: {e}")
    nested_combined_df = pd.concat(nested_dfs, ignore_index=True)
    print(f"Combined nested results have {len(nested_combined_df)} entries")

    # Step 2: Read all FASTA files from the genome_detective_fasta_folder
    fasta_files = glob.glob(os.path.join(genome_detective_fasta_folder, '*.fasta'))
    genome_detective_sequences = {}
    for fasta_file in fasta_files:
        sequences = read_fasta_sequences(fasta_file)
        genome_detective_sequences.update(sequences)
    print(f"Total genome detective sequences loaded: {len(genome_detective_sequences)}")

    # Step 3: Match names between CSV and FASTA to get combined data
    matched_names = []
    unmatched_sequences = []
    for index, row in nested_combined_df.iterrows():
        genome_detective_name = row['name'].strip()
        sequence = genome_detective_sequences.get(genome_detective_name, None)
        if sequence:
            matched_names.append((genome_detective_name, row['species'], row['assignment'], sequence))
        else:
            unmatched_sequences.append(genome_detective_name)
            print(f"Unmatched sequence for: {genome_detective_name}")

    print(f"Total matched names: {len(matched_names)}")
    print(f"Total unmatched names: {len(unmatched_sequences)}")

    # Step 4: Read unique passed sequences
    unique_passed_sequences = read_fasta_sequences(unique_passed_fasta_file)
    print(f"Total unique passed sequences loaded: {len(unique_passed_sequences)}")

    # Step 5: Match sequences and replace names with those from unique passed sequences
    final_data = []
    for genome_detective_name, species, assignment_column, sequence in matched_names:
        matched_name = next((name for name, seq in unique_passed_sequences.items() if seq == sequence), 'Unknown')
        final_data.append((matched_name, species, assignment_column, sequence))
        if matched_name == 'Unknown':
            print(f"Unmatched sequence for: {genome_detective_name}")

    # Step 6: Read BLAST genotyping results
    try:
        blast_df = pd.read_csv(blast_results_file, delimiter='\t', header=None, engine='python')
        blast_df.columns = ['name', 'taxonomy', 'score']
        blast_df['blast_serotype'] = blast_df['taxonomy'].apply(extract_serotype)
        blast_df['blast_genotype'] = blast_df['taxonomy'].apply(lambda x: x.split(';')[6] if len(x.split(';')) > 6 else 'Unknown')
        print(f"Loaded BLAST results with {len(blast_df)} entries")
    except Exception as e:
        print(f"Failed to read BLAST results file: {e}")
        sys.exit(1)

    # Step 7: Combine results
    final_combined_data = []
    for name, species, assignment_column, sequence in final_data:
        blast_entry = blast_df.loc[blast_df['name'] == name]
        blast_serotype = blast_entry['blast_serotype'].values[0] if not blast_entry.empty else 'Unknown'
        blast_genotype = blast_entry['blast_genotype'].values[0] if not blast_entry.empty else 'Unknown'
        genome_detective_serotype = extract_serotype(species)
        genome_detective_genotype = map_genotype(assignment_column)
        print(f"Name: {name}, GD Serotype: {genome_detective_serotype}, BLAST Serotype: {blast_serotype}, GD Genotype: {genome_detective_genotype}, BLAST Genotype: {blast_genotype}")
        final_combined_data.append((name, genome_detective_serotype, blast_serotype, genome_detective_genotype, blast_genotype, sequence))

    combined_df = pd.DataFrame(final_combined_data, columns=['name', 'genome_detective_serotype', 'blast_serotype', 'genome_detective_genotype', 'blast_genotype', 'sequence'])
    combined_df.to_csv(combined_output_file, index=False)
    print(f"Combined DataFrame has {len(combined_df)} matching entries")

    if len(combined_df) == 0:
        print("No matching entries found. Exiting.")
        return

    # Step 8: Generate confusion matrices
    def plot_confusion_matrix(true_labels, predicted_labels, unique_labels, title, output_file):
        cm = confusion_matrix(true_labels, predicted_labels, labels=unique_labels)
        report = classification_report(true_labels, predicted_labels, labels=unique_labels, target_names=unique_labels)
        print(report)

        report_file = os.path.join(os.path.dirname(output_file), f"classification_report_{title}.txt")
        with open(report_file, 'w') as f:
            f.write(report)

        plt.figure(figsize=(14, 10))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                    xticklabels=unique_labels, yticklabels=unique_labels,
                    cbar_kws={'shrink': 0.75}, linewidths=.5)
        plt.xlabel('Predicted Labels', fontsize=12)
        plt.ylabel('True Labels', fontsize=12)
        plt.title(f'Confusion Matrix - {title}', fontsize=15)
        plt.xticks(rotation=45, ha='right', fontsize=10)
        plt.yticks(rotation=0, fontsize=10)
        plt.tight_layout()
        plt.savefig(output_file)

    true_serotypes = combined_df['genome_detective_serotype']
    predicted_serotypes = combined_df['blast_serotype']
    unique_serotypes = sorted(true_serotypes.unique())

    if unique_serotypes:
        plot_confusion_matrix(true_serotypes, predicted_serotypes, unique_serotypes, 'Serotypes', confusion_matrix_output_file)

    true_genotypes = combined_df['genome_detective_genotype']
    predicted_genotypes = combined_df['blast_genotype']
    unique_genotypes = sorted(true_genotypes.unique())

    if unique_genotypes:
        genotype_confusion_matrix_output_file = confusion_matrix_output_file.replace('.png', '_genotypes.png')
        plot_confusion_matrix(true_genotypes, predicted_genotypes, unique_genotypes, 'Genotypes', genotype_confusion_matrix_output_file)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python generate_confusion_matrix.py <nested_results_folder> <genome_detective_fasta_folder> <blast_results_file> <unique_passed_fasta_file> <combined_output_file> <confusion_matrix_output_file>")
        sys.exit(1)
    
    nested_results_folder = sys.argv[1]
    genome_detective_fasta_folder = sys.argv[2]
    blast_results_file = sys.argv[3]
    unique_passed_fasta_file = sys.argv[4]
    combined_output_file = sys.argv[5]
    confusion_matrix_output_file = sys.argv[6]

    generate_confusion_matrix(nested_results_folder, genome_detective_fasta_folder, blast_results_file, unique_passed_fasta_file, combined_output_file, confusion_matrix_output_file)
