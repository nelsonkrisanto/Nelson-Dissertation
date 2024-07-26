import os
import sys
import pandas as pd
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score, balanced_accuracy_score, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import label_binarize
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
    print(f"Confusion Matrix: \n{cm}")

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
    primer_pairs = df['Primer Pairs'].unique()
    metrics = []

    all_labels = ['dengue_virus_type_1', 'dengue_virus_type_2', 'dengue_virus_type_3', 'dengue_virus_type_4', 'other_virus']

    for primer_pair in primer_pairs:
        subset = df[df['Primer Pairs'] == primer_pair]
        y_true = subset['Whole Genome Genotyping']
        y_pred = subset['Amplicons Genotyping']
        
        print(f"Calculating metrics for Primer Pair: {primer_pair}")
        print(f"True Labels: {y_true.values}")
        print(f"Predicted Labels: {y_pred.values}")
        
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
            'roc_auc': roc_auc
        })

    return pd.DataFrame(metrics)

def find_best_primer_pair(metrics_df):
    return metrics_df.sort_values(by='f1_score', ascending=False).iloc[0]

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python confusion_matrix_whole.py <path_to_csv_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    output_dir = 'results'  # Update this with the actual output directory

    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' does not exist.")
        sys.exit(1)

    df = pd.read_csv(file_path)
    metrics_df = process_primer_pairs(df, output_dir)
    
    best_primer_pair = find_best_primer_pair(metrics_df)
    
    output_file = f'{output_dir}/primer_pair_metrics.csv'
    metrics_df.to_csv(output_file, index=False)
    
    print("Best Primer Pair based on F1 score:")
    print(best_primer_pair)
    print(f"Metrics CSV saved to: {output_file}")
