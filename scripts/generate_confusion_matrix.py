import pandas as pd
import glob
from Bio import SeqIO
from sklearn.metrics import confusion_matrix, f1_score, precision_score, recall_score, classification_report
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

def process_species(species):
    if isinstance(species, str):
        return species.replace('dengue virus type ', 'DENV')
    return species

def process_assignment(assignment):
    if isinstance(assignment, str):
        assignment = assignment.split('_')[0]
        if assignment in ["dengue virus type 1", "dengue virus type 2", "dengue virus type 3", "dengue virus type 4"]:
            return "Unassigned"
        return assignment
    return assignment

def load_amplicon_data(path_pattern):
    all_files = glob.glob(path_pattern)
    df_list = []
    for file in all_files:
        df = pd.read_csv(file, sep=',')
        df_list.append(df)
    df_concat = pd.concat(df_list, ignore_index=True)
    return df_concat

def load_fasta_data(path_pattern):
    all_files = glob.glob(path_pattern)
    fasta_sequences = []
    for file in all_files:
        fasta_sequences.extend(list(SeqIO.parse(file, "fasta")))
    return fasta_sequences

def generate_confusion_matrix():
    amplicons_file_pattern_csv = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\GD\p*.csv'
    amplicons_file_pattern_fasta = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\GD\p*.fasta'
    clustered_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\test2\test2.csv'
    primer_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\filtered_linked_amplicons.csv'
    output_csv = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\output.csv'

    df_amplicons = load_amplicon_data(amplicons_file_pattern_csv)
    df_amplicons.columns = df_amplicons.columns.str.strip().str.replace('"', '')
    print("Amplicons Data Columns:", df_amplicons.columns)
    
    fasta_sequences = load_fasta_data(amplicons_file_pattern_fasta)
    
    fasta_df = pd.DataFrame([{"name": record.id, "sequence": str(record.seq)} for record in fasta_sequences])
    
    df_amplicons = df_amplicons.merge(fasta_df, on="name")
    
    df_clustered = pd.read_csv(clustered_file, sep=',')
    df_clustered.columns = df_clustered.columns.str.strip().str.replace('"', '')
    print("Clustered Data Columns:", df_clustered.columns)
    
    df_primers = pd.read_csv(primer_file)
    print("Primers Data Columns:", df_primers.columns)
    
    df_amplicons['species_amplicon'] = df_amplicons['species'].apply(process_species).astype(str)
    df_amplicons['assignment_amplicon'] = df_amplicons['assignment'].apply(process_assignment).astype(str)

    df_clustered['species_clustered'] = df_clustered['species'].apply(process_species).astype(str)
    df_clustered['assignment_clustered'] = df_clustered['assignment'].apply(process_assignment).astype(str)

    df_merged = df_primers.merge(df_amplicons, left_on='Sequence', right_on='sequence')
    
    df_merged = df_merged.drop(columns=['assignment', 'species', 'sequence'])

    df_merged = df_merged.merge(df_clustered, left_on='Accession_Number', right_on='name', suffixes=('_amplicon', '_clustered'))

    final_df = df_merged[['Amplicon_ID', 'Primer_Pair', 'Accession_Number', 'Sequence', 'species_clustered', 'species_amplicon', 'assignment_clustered', 'assignment_amplicon']]

    final_df['assignment_amplicon'].replace({"dengue virus type 1": "Unassigned", 
                                             "dengue virus type 2": "Unassigned", 
                                             "dengue virus type 3": "Unassigned", 
                                             "dengue virus type 4": "Unassigned"}, inplace=True)
    final_df['assignment_amplicon'].fillna("Unassigned", inplace=True)
    final_df['assignment_clustered'].fillna("Unassigned", inplace=True)

    final_df['species_amplicon'].replace("nan", "Unassigned", inplace=True)
    final_df['species_clustered'].replace("nan", "Unassigned", inplace=True)

    final_df.to_csv(output_csv, index=False)

    serotype_cm = confusion_matrix(final_df['species_clustered'], final_df['species_amplicon'])
    genotype_cm = confusion_matrix(final_df['assignment_clustered'], final_df['assignment_amplicon'])

    serotype_labels = sorted(final_df['species_clustered'].unique())
    genotype_labels = sorted(final_df['assignment_clustered'].unique())

    plt.figure(figsize=(12, 10))
    sns.heatmap(serotype_cm, annot=True, fmt='d', cmap='Blues', xticklabels=serotype_labels, yticklabels=serotype_labels)
    plt.title('Serotype Confusion Matrix')
    plt.xlabel('Amplicon Serotype')
    plt.ylabel('Clustered Serotype')
    plt.tight_layout()
    plt.savefig(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\serotype_confusion_matrix.png')
    plt.close()

    plt.figure(figsize=(20, 10))
    sns.heatmap(genotype_cm, annot=True, fmt='d', cmap='Blues', xticklabels=genotype_labels, yticklabels=genotype_labels)
    plt.title('Genotype Confusion Matrix')
    plt.xlabel('Amplicon Genotype')
    plt.ylabel('Clustered Genotype')
    plt.tight_layout()
    plt.savefig(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\genotype_confusion_matrix.png')
    plt.close()

    def calculate_metrics(df):
        y_true = df['species_clustered']
        y_pred = df['species_amplicon']
        f1 = f1_score(y_true, y_pred, average='weighted')
        precision = precision_score(y_true, y_pred, average='weighted')
        recall = recall_score(y_true, y_pred, average='weighted')
        report = classification_report(y_true, y_pred, output_dict=True)
        support = report['weighted avg']['support']
        return pd.Series({'F1_Score': f1, 'Precision': precision, 'Recall': recall, 'Support': support})

    overall_metrics = final_df.groupby('Primer_Pair').apply(calculate_metrics).reset_index()
    overall_metrics = overall_metrics.sort_values(by='F1_Score', ascending=False)
    overall_metrics.to_csv(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\metrics_overall.csv', index=False)

    best_f1_score = overall_metrics.iloc[0]
    with open(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\best_f1_score.txt', 'w') as file:
        file.write(f"Primer Pair with the Best F1 Score: {best_f1_score['Primer_Pair']} - F1 Score: {best_f1_score['F1_Score']}")

    serotypes = final_df['species_clustered'].unique()
    for serotype in serotypes:
        df_serotype = final_df[final_df['species_clustered'] == serotype]
        serotype_metrics = df_serotype.groupby('Primer_Pair').apply(calculate_metrics).reset_index()
        serotype_metrics = serotype_metrics.sort_values(by='F1_Score', ascending=False)
        serotype_metrics.to_csv(f'C:\\Users\\Nelso\\OneDrive\\Documents\\Thesis\\data\\conv\\metrics_serotype_{serotype}.csv', index=False)

    genotypes = final_df['assignment_clustered'].unique()
    for genotype in genotypes:
        df_genotype = final_df[final_df['assignment_clustered'] == genotype]
        genotype_metrics = df_genotype.groupby('Primer_Pair').apply(calculate_metrics).reset_index()
        genotype_metrics = genotype_metrics.sort_values(by='F1_Score', ascending=False)
        genotype_metrics.to_csv(f'C:\\Users\\Nelso\\OneDrive\\Documents\\Thesis\\data\\conv\\metrics_genotype_{genotype}.csv', index=False)

    # Format the mismatched classifications output
    mismatched_classifications = final_df[final_df['assignment_clustered'] != final_df['assignment_amplicon']]
    mismatched_classifications = mismatched_classifications.groupby(['Primer_Pair', 'assignment_clustered', 'assignment_amplicon']).size().reset_index(name='Number of observations')
    mismatched_classifications.columns = ['Primer Pair', 'Clustered Genotype', 'Amplicon Genotype', 'Number of Observations']
    mismatched_classifications = mismatched_classifications.sort_values(by='Number of Observations', ascending=False)
    mismatched_classifications.to_csv(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\mismatched_classifications.csv', index=False)

    # Format the missed classifications output
    missed_classifications = final_df[final_df['assignment_clustered'] != final_df['assignment_amplicon']]
    missed_classifications = missed_classifications.groupby(['Primer_Pair', 'assignment_clustered']).size().reset_index(name='Number of observations')
    missed_classifications.columns = ['Primer Pair', 'Clustered Genotype', 'Number of Observations']
    missed_classifications = missed_classifications.sort_values(by='Number of Observations', ascending=False)
    missed_classifications.to_csv(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\missed_classifications.csv', index=False)

    # Statistical test for misclassification
    misclassified_genotypes = missed_classifications['Clustered Genotype'].value_counts().reset_index()
    misclassified_genotypes.columns = ['Genotype', 'Count']
    chi2_missed, p_missed, _, _ = chi2_contingency([misclassified_genotypes['Count'], [misclassified_genotypes['Count'].sum() - x for x in misclassified_genotypes['Count']]])
    
    # Statistical test for mismatched classification
    mismatched_genotypes = mismatched_classifications['Clustered Genotype'].value_counts().reset_index()
    mismatched_genotypes.columns = ['Genotype', 'Count']
    chi2_mismatched, p_mismatched, _, _ = chi2_contingency([mismatched_genotypes['Count'], [mismatched_genotypes['Count'].sum() - x for x in mismatched_genotypes['Count']]])
    
    with open(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\statistical_test_results.txt', 'w') as file:
        file.write(f"Chi-square Test for Missed Genotypes\n")
        file.write(f"Chi-square value: {chi2_missed}\n")
        file.write(f"P-value: {p_missed}\n")
        if p_missed < 0.05:
            file.write("There is a statistically significant difference in the misclassification of genotypes.\n")
        else:
            file.write("There is no statistically significant difference in the misclassification of genotypes.\n")
        
        file.write("\nChi-square Test for Mismatched Genotypes\n")
        file.write(f"Chi-square value: {chi2_mismatched}\n")
        file.write(f"P-value: {p_mismatched}\n")
        if p_mismatched < 0.05:
            file.write("There is a statistically significant difference in the mismatched classification of genotypes.\n")
        else:
            file.write("There is no statistically significant difference in the mismatched classification of genotypes.\n")

if __name__ == "__main__":
    generate_confusion_matrix()
