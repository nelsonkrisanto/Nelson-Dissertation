import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
missed_classifications = pd.read_csv(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\missed_classifications.csv')
mismatched_classifications = pd.read_csv(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\mismatched_classifications.csv')

# Combine all relevant data into a single DataFrame for mismatched classifications
mismatched_data = []

for idx, row in mismatched_classifications.iterrows():
    mismatched_data.append({
        'Type': 'Mismatched',
        'Primer Pair': row['Primer Pair'],
        'Clustered': row['Clustered Genotype'],
        'Amplicon': row['Amplicon Genotype'],
        'Number of Observations': row['Number of Observations']
    })

mismatched_df = pd.DataFrame(mismatched_data)

# Prepare data for mismatched heatmap
mismatched_heatmap_data = mismatched_df.pivot_table(
    index=['Primer Pair', 'Clustered'],
    columns='Amplicon',
    values='Number of Observations',
    aggfunc='sum',
    fill_value=0
)

# Plot mismatched heatmap
plt.figure(figsize=(20, 15))
sns.heatmap(mismatched_heatmap_data, annot=True, fmt='d', cmap='rocket', linewidths=.5)
plt.title('Heatmap of Mismatched Classifications')
plt.xlabel('Amplicon Genotype')
plt.ylabel('Primer Pair - Clustered Genotype')
plt.tight_layout()
plt.savefig(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\mismatched_heatmap.png')
plt.show()

# Combine all relevant data into a single DataFrame for missed classifications
missed_data = []

for idx, row in missed_classifications.iterrows():
    missed_data.append({
        'Type': 'Missed',
        'Primer Pair': row['Primer Pair'],
        'Clustered': row['Clustered Genotype'],
        'Genotype': row['Clustered Genotype'],
        'Number of Observations': row['Number of Observations']
    })

missed_df = pd.DataFrame(missed_data)

# Prepare data for missed heatmap
missed_heatmap_data = missed_df.pivot_table(
    index=['Primer Pair', 'Clustered'],
    columns='Genotype',
    values='Number of Observations',
    aggfunc='sum',
    fill_value=0
)

# Plot missed heatmap
plt.figure(figsize=(20, 15))
sns.heatmap(missed_heatmap_data, annot=True, fmt='d', cmap='rocket', linewidths=.5)
plt.title('Heatmap of Missed Classifications')
plt.xlabel('Missed Genotype')
plt.ylabel('Primer Pair - Clustered Genotype')
plt.tight_layout()
plt.savefig(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\conv\missed_heatmap.png')
plt.show()

# Print to verify
print("Mismatched Heatmap Data:\n", mismatched_heatmap_data)
print("Missed Heatmap Data:\n", missed_heatmap_data)
