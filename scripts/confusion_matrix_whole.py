import os
import pandas as pd
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score, balanced_accuracy_score, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import label_binarize
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def plot_confusion_matrix(cm, labels, title, output_file):
    plt.figure(figsize=(10, 7))
    sns.heatmap(cm, annot=True, fmt='d', xticklabels=labels, yticklabels=labels, cmap='Blues')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title(title)
    plt.savefig(output_file)
    plt.close()

def calculate_metrics(y_true, y_pred, labels):
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    precision = precision_score(y_true, y_pred, average='macro', zero_division=0)
    recall = recall_score(y_true, y_pred, average='macro', zero_division=0)
    f1 = f1_score(y_true, y_pred, average='macro', zero_division=0)
    specificity = []
    npv = []
    for i, label in enumerate(labels):
        tn = cm.sum() - (cm[i, :].sum() + cm[:, i].sum() - cm[i, i])
        fp = cm[:, i].sum() - cm[i, i]
        fn = cm[i, :].sum() - cm[i, i]
        tn_fp = tn + fp
        tn_fn = tn + fn
        specificity.append(tn / tn_fp if tn_fp != 0 else 0)
        npv.append(tn / tn_fn if tn_fn != 0 else 0)
    specificity = np.mean(specificity)
    npv = np.mean(npv)
    mcc = matthews_corrcoef(y_true, y_pred)
    balanced_acc = balanced_accuracy_score(y_true, y_pred)
    
    y_true_bin = label_binarize(y_true, classes=labels)
    y_pred_bin = label_binarize(y_pred, classes=labels)
    try:
        roc_auc = roc_auc_score(y_true_bin, y_pred_bin, average='macro')
    except ValueError:
        roc_auc = 'N/A'

    return precision, recall, f1, specificity, npv, mcc, balanced_acc, roc_auc, cm

def process_primer_pairs(df, output_dir):
    primer_pairs = df['Primer_Pair'].unique()
    metrics = []

    all_labels = ['dengue_virus_type_1', 'dengue_virus_type_2', 'dengue_virus_type_3', 'dengue_virus_type_4', 'other_virus']

    for primer_pair in primer_pairs:
        subset = df[df['Primer_Pair'] == primer_pair]
        y_true = subset['Whole_Genome_Genotyping']
        y_pred = subset['Amplicons_Genotyping']
        
        precision, recall, f1, specificity, npv, mcc, balanced_acc, roc_auc, cm = calculate_metrics(y_true, y_pred, all_labels)
        
        plot_confusion_matrix(cm, all_labels, f'Confusion Matrix for {primer_pair}', f'{output_dir}/{primer_pair}_confusion_matrix.png')
        
        metrics.append({
            'primer_pair': primer_pair,
            'precision': precision,
            'recall': recall,
            'f1_score': f1,
            'specificity': specificity,
            'npv': npv,
            'mcc': mcc,
            'balanced_accuracy': balanced_acc,
            'roc_auc': roc_auc,
            'amplicon_length': subset['Amplicon_Length'].mean()  # assuming 'Amplicon_Length' column exists
        })

    return pd.DataFrame(metrics)

def find_best_primer_pair(metrics_df):
    return metrics_df.sort_values(by='f1_score', ascending=False).iloc[0]

def breakdown_by_serotype(df, serotype, output_dir):
    subset = df[df['Whole_Genome_Genotyping'] == serotype]
    primer_pairs = subset['Primer_Pair'].unique()
    metrics = []

    all_labels = ['dengue_virus_type_1', 'dengue_virus_type_2', 'dengue_virus_type_3', 'dengue_virus_type_4', 'other_virus']

    for primer_pair in primer_pairs:
        subset_pair = subset[subset['Primer_Pair'] == primer_pair]
        y_true = subset_pair['Whole_Genome_Genotyping']
        y_pred = subset_pair['Amplicons_Genotyping']
        
        precision, recall, f1, specificity, npv, mcc, balanced_acc, roc_auc, cm = calculate_metrics(y_true, y_pred, all_labels)
        
        plot_confusion_matrix(cm, all_labels, f'Confusion Matrix for {primer_pair} ({serotype})', f'{output_dir}/{primer_pair}_{serotype}_confusion_matrix.png')
        
        metrics.append({
            'primer_pair': primer_pair,
            'serotype': serotype,
            'precision': precision,
            'recall': recall,
            'f1_score': f1,
            'specificity': specificity,
            'npv': npv,
            'mcc': mcc,
            'balanced_accuracy': balanced_acc,
            'roc_auc': roc_auc,
            'amplicon_length': subset_pair['Amplicon_Length'].mean()  # assuming 'Amplicon_Length' column exists
        })

    return pd.DataFrame(metrics)

def overall_confusion_matrix(df, output_dir):
    y_true = df['Whole_Genome_Genotyping']
    y_pred = df['Amplicons_Genotyping']
    
    all_labels = ['dengue_virus_type_1', 'dengue_virus_type_2', 'dengue_virus_type_3', 'dengue_virus_type_4', 'other_virus']
    
    precision, recall, f1, specificity, npv, mcc, balanced_acc, roc_auc, cm = calculate_metrics(y_true, y_pred, all_labels)
    
    plot_confusion_matrix(cm, all_labels, 'Overall Confusion Matrix for Dengue Virus', f'{output_dir}/overall_confusion_matrix.png')
    
    return pd.DataFrame([{
        'primer_pair': 'Overall',
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'specificity': specificity,
        'npv': npv,
        'mcc': mcc,
        'balanced_accuracy': balanced_acc,
        'roc_auc': roc_auc,
        'amplicon_length': df['Amplicon_Length'].mean()  # assuming 'Amplicon_Length' column exists
    }])

def correlation_analysis(metrics_df, output_dir):
    correlation = metrics_df[['f1_score', 'amplicon_length']].corr().iloc[0, 1]
    plt.figure(figsize=(10, 7))
    sns.regplot(x='amplicon_length', y='f1_score', data=metrics_df)
    plt.title(f'Correlation between F1-score and Amplicon Length: {correlation:.2f}')
    plt.xlabel('Amplicon Length')
    plt.ylabel('F1-score')
    plt.savefig(f'{output_dir}/correlation_f1score_amplicon_length.png')
    plt.close()

def commonly_misclassified(df, output_dir):
    misclassified = df[df['Whole_Genome_Genotyping'] != df['Amplicons_Genotyping']]
    misclassified_summary = misclassified.groupby(['Whole_Genome_Genotyping', 'Amplicons_Genotyping']).size().reset_index(name='count')
    misclassified_summary = misclassified_summary.sort_values(by='count', ascending=False)
    misclassified_summary.to_csv(f'{output_dir}/commonly_misclassified_genotypes.csv', index=False)
    print("Commonly misclassified genotypes:")
    print(misclassified_summary.head(10))  # print top 10 misclassified genotypes

if __name__ == "__main__":
    file_path = "C:/Users/Nelso/OneDrive/Documents/Thesis/data/conv/conv_whole_confusion_matrix.csv"
    output_dir = "C:/Users/Nelso/OneDrive/Documents/Thesis/results/conventional"

    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' does not exist.")
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df = pd.read_csv(file_path)

    # Check columns for debugging
    print("Columns in DataFrame:", df.columns.tolist())

    # Overall confusion matrix
    overall_metrics = overall_confusion_matrix(df, output_dir)
    
    # Process by primer pairs
    metrics_df = process_primer_pairs(df, output_dir)
    
    # Breakdown by serotype
    serotypes = ['dengue_virus_type_1', 'dengue_virus_type_2', 'dengue_virus_type_3', 'dengue_virus_type_4']
    for serotype in serotypes:
        serotype_metrics = breakdown_by_serotype(df, serotype, output_dir)
        metrics_df = pd.concat([metrics_df, serotype_metrics], ignore_index=True)
    
    # Correlation analysis
    correlation_analysis(metrics_df, output_dir)

    # Commonly misclassified genotypes
    commonly_misclassified(df, output_dir)

    # Save metrics
    output_file = f'{output_dir}/primer_pair_metrics.csv'
    metrics_df.to_csv(output_file, index=False)

    # Find the best primer pair based on F1 score
    best_primer_pair = find_best_primer_pair(metrics_df)

    print("Best Primer Pair based on F1 score:")
    print(best_primer_pair)
    print(f"Metrics CSV saved to: {output_file}")
