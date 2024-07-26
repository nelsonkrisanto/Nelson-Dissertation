import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Load the data from the TSV file
file_path = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\individual_conv_primer_details.tsv'
df = pd.read_csv(file_path, delimiter='\t')

# Display the dataframe to check the loaded data
print(df.head())
print("Columns in DataFrame:", df.columns)

# Define colors for the regions and genotypes
region_colors = {'NS1': '#EDD9BA', 'NS3': '#F0E2B6', 'NS5': '#DDDAF4'}
genotype_colors = {
    'DENV1': '#FF5733',  # Strong red
    'DENV2': '#33C1FF',  # Strong blue
    'DENV3': '#75FF33',  # Strong green
    'DENV4': '#FF33EC',  # Strong pink
    'ALL': '#000000'     # Black for 'ALL'
}

# Plot settings
fig, ax = plt.subplots(figsize=(12, 3))
regions = {'NS1': (1000, 2000), 'NS3': (3000, 4000), 'NS5': (5000, 7000)}

# Plot regions
for region, (start, end) in regions.items():
    ax.add_patch(patches.Rectangle((start, 0), end - start, 1, edgecolor='black', facecolor=region_colors[region], alpha=0.3, label=region))

# Group data by region and genotype, then calculate percentages
grouped = df.groupby(['Region', 'Genotype']).size().unstack(fill_value=0)
percentages = grouped.divide(grouped.sum(axis=1), axis=0) * 100

# Plot percentage and population number for each genotype in each region
for region, (start, end) in regions.items():
    region_data = percentages.loc[region]
    total = region_data.sum()
    y_position = 0.8  # Starting y position
    for genotype, percentage in region_data.items():
        count = grouped.loc[region, genotype]
        if count > 0:
            ax.text((start + end) / 2, y_position, f"• {genotype}: {count} ({percentage:.2f}%)", fontsize=8, color=genotype_colors[genotype], ha='center', va='center')
            y_position -= 0.15  # Adjust y position for the next line

# Set plot limits and labels
ax.set_xlim(0, 10707)
ax.set_ylim(0, 1)
ax.set_yticks([])

# Remove y-axis label
ax.set_ylabel('')

# Remove x-axis label and set custom labels
ax.set_xticks([1500, 3500, 6000])
ax.set_xticklabels(['NS1', 'NS3', 'NS5'])

# Create a custom legend for genotypes
handles = [patches.Patch(color=genotype_colors['ALL'], label='ALL'),
           patches.Patch(color=genotype_colors['DENV1'], label='DENV1'),
           patches.Patch(color=genotype_colors['DENV2'], label='DENV2'),
           patches.Patch(color=genotype_colors['DENV3'], label='DENV3'),
           patches.Patch(color=genotype_colors['DENV4'], label='DENV4')]

ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.1, 1))

plt.tight_layout()

# Save the plot as an image file
output_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\regions_heatmap_labeled_v4.png'
plt.savefig(output_file, format='png', dpi=300)

plt.show()
