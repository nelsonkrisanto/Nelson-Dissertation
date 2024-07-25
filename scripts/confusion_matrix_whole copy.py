# Use the Agg backend for non-GUI environments
import matplotlib
matplotlib.use('Agg')

import sys
import os
import pandas as pd
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import warnings
from sklearn.exceptions import UndefinedMetricWarning

# Suppress UndefinedMetricWarning
warnings.filterwarnings(action='ignore', category=UndefinedMetricWarning)

# Function to load data
def load_data(file_path):
    return pd.read_csv(file_path)

# Function to plot confusion matrix
def plot_confusion_matrix(cm, labels, title, output_file):
    plt.figure(figsize=(10, 7))
    sns.heatmap(cm, annot=True, fmt='d', xticklabels=labels, yticklabels=labels, cmap='Blues')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title(title)
    plt.savefig(output_file)
    plt.close()

# Function to calculate metrics
def calculate_metrics(df, output_dir):
    primer_pairs = df['Primer Pairs'].unique()
    metrics = []

    for primer_pair in primer_pairs:
        subset = df[df['Primer Pairs'] == primer_pair]
        y_true = subset['Whole Genome Genotyping']
        y_pred = subset['Amplicons Genotyping']
        
        # Calculate confusion matrix
        cm = confusion_matrix(y_true, y_pred, labels=y_true.unique())
        plot_confusion_matrix(cm, y_true.unique(), f'Confusion Matrix for {primer_pair}', os.path.join(output_dir, f'{primer_pair}_confusion_matrix.png'))
        
        # Calculate precision, recall, and F1 score, handle zero division manually
        precision = precision_score(y_true, y_pred, average='macro') if np.sum(cm.sum(axis=0)) > 0 else 0.0
        recall = recall_score(y_true, y_pred, average='macro') if np.sum(cm.sum(axis=1)) > 0 else 0.0
        f1 = f1_score(y_true, y_pred, average='macro') if precision > 0 and recall > 0 else 0.0
        
        metrics.append({
            'primer_pair': primer_pair,
            'precision': precision,
            'recall': recall,
            'f1_score': f1
        })

    return pd.DataFrame(metrics)

# Function to find the best primer pair
def find_best_primer_pair(metrics_df):
    return metrics_df.sort_values(by='f1_score', ascending=False).iloc[0]

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python confusion_matrix_whole.py <file_path>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    output_dir = os.path.dirname(file_path)
    
    # Load the dataset
    df = load_data(file_path)
    
    # Calculate metrics for each primer pair
    metrics_df = calculate_metrics(df, output_dir)
    
    # Find the best primer pair based on F1 score
    best_primer_pair = find_best_primer_pair(metrics_df)
    
    # Save the metrics to a CSV file
    output_file = os.path.join(output_dir, 'primer_pair_metrics.csv')
    metrics_df.to_csv(output_file, index=False)
    
    # Print the best primer pair
    print("Best Primer Pair based on F1 score:")
    print(best_primer_pair)
    print(f"Metrics CSV saved to: {output_file}")
