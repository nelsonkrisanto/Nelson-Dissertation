import os
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix, precision_recall_fscore_support, balanced_accuracy_score, matthews_corrcoef

def plot_confusion_matrix(cm, labels, title, output_file):
    plt.figure(figsize=(14, 12))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=labels, yticklabels=labels)
    plt.title(title)
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def calculate_metrics(y_true, y_pred, labels):
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    precision, recall, f1, _ = precision_recall_fscore_support(y_true, y_pred, average='weighted', labels=labels, zero_division=0)
    specificity = np.mean([cm[i, i] / cm[i, :].sum() if cm[i, :].sum() != 0 else 0 for i in range(len(labels))])
    npv = np.mean([cm[i, i] / cm[:, i].sum() if cm[:, i].sum() != 0 else 0 for i in range(len(labels))])
    mcc = matthews_corrcoef(y_true, y_pred)
    balanced_acc = balanced_accuracy_score(y_true, y_pred)
    roc_auc = 'N/A'
    return precision, recall, f1, specificity, npv, mcc, balanced_acc, roc_auc, cm

def standardize_labels(df):
    conversion_dict = {
        '1I': 'DENV-1_genotype_I', '1II': 'DENV-1_genotype_II', '1III': 'DENV-1_genotype_III', '1IV': 'DENV-1_genotype_IV', '1V': 'DENV-1_genotype_V',
        '2I': 'DENV-2_genotype_Asian_I', '2II': 'DENV-2_genotype_Asian_II', '2III': 'DENV-2_genotype_American', '2IV': 'DENV-2_genotype_Cosmopolitan',
        '3I': 'DENV-3_genotype_I', '3II': 'DENV-3_genotype_II', '3III': 'DENV-3_genotype_III', '3IV': 'DENV-3_genotype_IV', '3V': 'DENV-3_genotype_V',
        '4I': 'DENV-4_genotype_I', '4II': 'DENV-4_genotype_II'
    }

    # Standardize Genome Detective Genotype
    df['Genome Detective Genotype'] = df['Genome Detective Genotype'].apply(lambda x: conversion_dict.get(x.split('_')[0], 'Unassigned') if x.split('_')[0] in conversion_dict else 'Unassigned')

    # Standardize Genome Detective Serotype and Blastx Serotype
    df['Genome Detective Serotype'] = df['Genome Detective Serotype'].replace({'dengue virus type 1': 'DENV1', 'dengue virus type 2': 'DENV2', 'dengue virus type 3': 'DENV3', 'dengue virus type 4': 'DENV4', '0': 'Unassigned'})
    df['Blastx Serotype'] = df['Blastx Serotype'].replace({'dengue_virus_type_1': 'DENV1', 'dengue_virus_type_2': 'DENV2', 'dengue_virus_type_3': 'DENV3', 'dengue_virus_type_4': 'DENV4', 'undetermined': 'Unassigned'})

    # Standardize Blastx Genotype
    df['Blastx Genotype'] = df['Blastx Genotype'].replace({'undetermined': 'Unassigned'})

    return df

def process_data(df):
    df = standardize_labels(df)

    # Separate y_true and y_pred for Genotype and Serotype
    y_true_geno = df['Genome Detective Genotype']
    y_pred_geno = df['Blastx Genotype']

    y_true_sero = df['Genome Detective Serotype']
    y_pred_sero = df['Blastx Serotype']

    return y_true_geno, y_pred_geno, y_true_sero, y_pred_sero

def save_metrics_to_txt(file_path, metrics_dict):
    with open(file_path, 'w') as file:
        for key, value in metrics_dict.items():
            file.write(f"{key}: {value}\n")

if __name__ == "__main__":
    file_path = r"C:\Users\Nelso\OneDrive\Documents\Thesis\data\test2\VS65.csv"
    output_dir = r"C:\Users\Nelso\OneDrive\Documents\Thesis\results\test2\65"

    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' does not exist.")
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df = pd.read_csv(file_path)

    y_true_geno, y_pred_geno, y_true_sero, y_pred_sero = process_data(df)

    all_labels_geno = np.union1d(y_true_geno.unique(), y_pred_geno.unique())
    all_labels_sero = np.union1d(y_true_sero.unique(), y_pred_sero.unique())

    precision_geno, recall_geno, f1_geno, specificity_geno, npv_geno, mcc_geno, balanced_acc_geno, roc_auc_geno, cm_geno = calculate_metrics(y_true_geno, y_pred_geno, all_labels_geno)
    precision_sero, recall_sero, f1_sero, specificity_sero, npv_sero, mcc_sero, balanced_acc_sero, roc_auc_sero, cm_sero = calculate_metrics(y_true_sero, y_pred_sero, all_labels_sero)

    print("Genotype Metrics:")
    print(f"Precision: {precision_geno}")
    print(f"Recall: {recall_geno}")
    print(f"F1 Score: {f1_geno}")
    print(f"Specificity: {specificity_geno}")
    print(f"NPV: {npv_geno}")
    print(f"MCC: {mcc_geno}")
    print(f"Balanced Accuracy: {balanced_acc_geno}")
    print(f"ROC AUC: {roc_auc_geno}")

    print("Serotype Metrics:")
    print(f"Precision: {precision_sero}")
    print(f"Recall: {recall_sero}")
    print(f"F1 Score: {f1_sero}")
    print(f"Specificity: {specificity_sero}")
    print(f"NPV: {npv_sero}")
    print(f"MCC: {mcc_sero}")
    print(f"Balanced Accuracy: {balanced_acc_sero}")
    print(f"ROC AUC: {roc_auc_sero}")

    # Save metrics to text file
    metrics_file_path = os.path.join(output_dir, 'metrics.txt')
    metrics_dict = {
        'Genotype Precision': precision_geno,
        'Genotype Recall': recall_geno,
        'Genotype F1 Score': f1_geno,
        'Genotype Specificity': specificity_geno,
        'Genotype NPV': npv_geno,
        'Genotype MCC': mcc_geno,
        'Genotype Balanced Accuracy': balanced_acc_geno,
        'Genotype ROC AUC': roc_auc_geno,
        'Serotype Precision': precision_sero,
        'Serotype Recall': recall_sero,
        'Serotype F1 Score': f1_sero,
        'Serotype Specificity': specificity_sero,
        'Serotype NPV': npv_sero,
        'Serotype MCC': mcc_sero,
        'Serotype Balanced Accuracy': balanced_acc_sero,
        'Serotype ROC AUC': roc_auc_sero
    }
    save_metrics_to_txt(metrics_file_path, metrics_dict)

    plot_confusion_matrix(cm_geno, all_labels_geno, 'Confusion Matrix for Genotypes', f'{output_dir}/confusion_matrix_genotypes.png')
    plot_confusion_matrix(cm_sero, all_labels_sero, 'Confusion Matrix for Serotypes', f'{output_dir}/confusion_matrix_serotypes.png')
