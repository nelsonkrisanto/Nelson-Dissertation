import sys

# Use the Agg backend for non-GUI environments
import matplotlib
matplotlib.use('Agg')

import pandas as pd
from sklearn.metrics import confusion_matrix, classification_report
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
file_path = sys.argv[1]
print(f"Reading data from {file_path}")
try:
    df = pd.read_csv(file_path)
    print("CSV Data loaded successfully")
except Exception as e:
    print(f"Failed to read CSV file: {e}")
    sys.exit(1)

print("CSV Data:")
print(df)

# Clean the data
df = df.dropna(subset=['Genome Detective', 'Genotype Assigned Blastx'])
df['Genome Detective'] = df['Genome Detective'].astype(str)
df['Genotype Assigned Blastx'] = df['Genotype Assigned Blastx'].astype(str)

print("Cleaned CSV Data:")
print(df)

# Extract the true labels and predicted labels
try:
    true_labels = df['Genome Detective']
    predicted_labels = df['Genotype Assigned Blastx']
    print("True Labels extracted:")
    print(true_labels)
    print("Predicted Labels extracted:")
    print(predicted_labels)
except KeyError as e:
    print(f"Missing expected column in CSV file: {e}")
    sys.exit(1)

# Create a confusion matrix
try:
    cm = confusion_matrix(true_labels, predicted_labels, labels=true_labels.unique())
    print("Confusion Matrix created:")
    print(cm)
except Exception as e:
    print(f"Failed to create confusion matrix: {e}")
    sys.exit(1)

# Print classification report
print("\nClassification Report:")
report = classification_report(true_labels, predicted_labels, labels=true_labels.unique(), target_names=true_labels.unique())
print(report)

# Save the classification report to a file
output_dir = "/home/people/23203786/scratch/Nelson-Dissertation/results/confusion_matrix/"
report_file = output_dir + "classification_report.txt"
try:
    with open(report_file, 'w') as f:
        f.write(report)
    print(f"Classification report saved to {report_file}")
except Exception as e:
    print(f"Failed to save classification report: {e}")
    sys.exit(1)

# Display the confusion matrix using seaborn
plt.figure(figsize=(14, 10))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
            xticklabels=true_labels.unique(), yticklabels=true_labels.unique(),
            cbar_kws={'shrink': 0.75}, linewidths=.5)
plt.xlabel('Predicted Labels', fontsize=12)
plt.ylabel('True Labels', fontsize=12)
plt.title('Confusion Matrix', fontsize=15)
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.yticks(rotation=0, fontsize=10)
plt.tight_layout()

# Save the confusion matrix as an image file
output_file = output_dir + "confusion_matrix.png"
try:
    plt.savefig(output_file)
    print(f"Confusion matrix saved to {output_file}")
except Exception as e:
    print(f"Failed to save confusion matrix image: {e}")
    sys.exit(1)

# Optionally show the plot
# plt.show()
