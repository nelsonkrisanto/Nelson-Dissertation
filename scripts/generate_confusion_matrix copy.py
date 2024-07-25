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
            matched_names.append((genome_detective_name, row['species'], sequence))
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
    for genome_detective_name, species, sequence in matched_names:
        matched_name = next((name for name, seq in unique_passed_sequences.items() if seq == sequence), 'Unknown')
        final_data.append((matched_name, species, sequence))
        if matched_name == 'Unknown':
            print(f"Unmatched sequence for: {genome_detective_name}")

    # Step 6: Read BLAST genotyping results
    try:
        blast_df = pd.read_csv(blast_results_file, delimiter='\t', header=None, engine='python')
        blast_df.columns = ['name', 'taxonomy', 'score']
        blast_df['blast_serotype'] = blast_df['taxonomy'].apply(extract_serotype)
        print(f"Loaded BLAST results with {len(blast_df)} entries")
    except Exception as e:
        print(f"Failed to read BLAST results file: {e}")
        sys.exit(1)

    # Step 7: Combine results
    final_combined_data = []
    for name, species, sequence in final_data:
        blast_serotype = blast_df.loc[blast_df['name'] == name, 'blast_serotype'].values
        blast_serotype = blast_serotype[0] if len(blast_serotype) > 0 else 'Unknown'
        genome_detective_serotype = extract_serotype(species)
        final_combined_data.append((name, genome_detective_serotype, blast_serotype, sequence))

    combined_df = pd.DataFrame(final_combined_data, columns=['name', 'genome_detective_serotype', 'blast_serotype', 'sequence'])
    combined_df.to_csv(combined_output_file, index=False)
    print(f"Combined DataFrame has {len(combined_df)} matching entries")

    if len(combined_df) == 0:
        print("No matching entries found. Exiting.")
        return

    # Step 8: Generate confusion matrix
    true_labels = combined_df['genome_detective_serotype']
    predicted_labels = combined_df['blast_serotype']
    unique_labels = sorted(true_labels.unique())

    if not unique_labels:
        print("No unique labels found. Exiting.")
        return

    cm = confusion_matrix(true_labels, predicted_labels, labels=unique_labels)

    report = classification_report(true_labels, predicted_labels, labels=unique_labels, target_names=unique_labels)
    print(report)

    report_file = os.path.join(os.path.dirname(confusion_matrix_output_file), "classification_report.txt")
    with open(report_file, 'w') as f:
        f.write(report)

    plt.figure(figsize=(14, 10))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=unique_labels, yticklabels=unique_labels,
                cbar_kws={'shrink': 0.75}, linewidths=.5)
    plt.xlabel('Predicted Labels', fontsize=12)
    plt.ylabel('True Labels', fontsize=12)
    plt.title('Confusion Matrix', fontsize=15)
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.yticks(rotation=0, fontsize=10)
    plt.tight_layout()

    plt.savefig(confusion_matrix_output_file)

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
