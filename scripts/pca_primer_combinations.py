import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

# Function to draw confidence ellipse
def confidence_ellipse(x, y, ax, n_std=1.96, facecolor='none', **kwargs):
    if x.size != y.size:
        raise ValueError("x and y must be the same size")
    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, facecolor=facecolor, **kwargs)
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)
    transf = (transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y))
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

# Load the input data
input_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\individual_conv_primer_combination_details.tsv'
df = pd.read_csv(input_file, delimiter='\t')

# Display the dataframe to check the loaded data
print(df.head())
print("Columns in DataFrame:", df.columns)

# Store the relevant information for later use
info_columns = ['Combination_Name', 'Genotype_Forward', 'Genotype_Reverse', 'Forward_Primer', 'Reverse_Primer']
info_df = df[info_columns + ['Region']].copy()

# Rename the columns
info_df.rename(columns={'Genotype_Forward': 'Serotype_Forward', 'Genotype_Reverse': 'Serotype_Reverse'}, inplace=True)

# Concatenate columns to create unique row keys
df['Row_Key'] = df.apply(lambda row: f"{row['Combination_Name']}_{row['Region']}", axis=1)
info_df['Row_Key'] = df['Row_Key']  # Add Row_Key to info_df for later merging
df.set_index('Row_Key', inplace=True)

# Drop the columns with leading information and any other unnecessary columns
df.drop(columns=info_columns + ['Reference'], inplace=True)

# Remove columns that start with 'Sequence_'
df = df[[col for col in df.columns if not col.startswith('Sequence_')]]

# Handle missing values by filling them only for numeric columns
numeric_columns = df.select_dtypes(include=['number']).columns
df[numeric_columns] = df[numeric_columns].fillna(df[numeric_columns].mean())

# Convert categorical columns to category dtype and then to dummy variables without dropping the first category
categorical_columns = df.select_dtypes(include=['object']).columns
df = pd.get_dummies(df, columns=categorical_columns, drop_first=False)

# Verify if the region columns are included
print("Columns after get_dummies:\n", df.columns)

# Standardize the data
scaler = StandardScaler()
df_scaled = scaler.fit_transform(df)

# Perform PCA
pca = PCA(n_components=2)
principal_components = pca.fit_transform(df_scaled)

# Calculate explained variance ratio
explained_variance_ratio = pca.explained_variance_ratio_

# Manually calculate cos2 for each observation
row_coords_squared = pd.DataFrame(principal_components ** 2, index=df.index)
row_cos2 = row_coords_squared.div(row_coords_squared.sum(axis=1), axis=0)
row_cos2.columns = ['cos2_dim1', 'cos2_dim2']

# Manually calculate contrib for each variable
loadings = pca.components_.T
eigenvalues = pca.explained_variance_
squared_loadings = loadings ** 2
contrib = squared_loadings * eigenvalues[np.newaxis, :]
contrib = pd.DataFrame(contrib.T, index=['contrib_dim1', 'contrib_dim2'], columns=df.columns).T

# Create a DataFrame for the principal components
famd_row_coords = pd.DataFrame(principal_components, columns=['Dim1', 'Dim2'], index=df.index)

# Merge the coordinates with the serotype and primer name information
famd_row_coords = famd_row_coords.merge(info_df, left_index=True, right_on='Row_Key')

# Enhance the plot
plt.figure(figsize=(14, 10))

# Define custom colors for each serotype
color_map = {
    'DENV1': '#FF5733',  # Strong red
    'DENV2': '#33C1FF',  # Strong blue
    'DENV3': '#75FF33',  # Strong green
    'DENV4': '#FF33EC',  # Strong pink
    'ALL': '#000000'     # Black for 'ALL'
}

# Define custom colors for each region
region_color_map = {
    'NS1': '#EDD9BA',
    'NS3': '#F0E2B6',
    'NS5': '#DDDAF4'
}

# Use the 'Serotype_Forward' column for classification (or 'Serotype_Reverse' if preferred)
serotypes = famd_row_coords['Serotype_Forward']
unique_serotypes = serotypes.unique()

# Scatter plot with colors
ax = plt.gca()
for serotype in unique_serotypes:
    idx = serotypes == serotype
    plt.scatter(famd_row_coords.loc[idx, 'Dim1'], famd_row_coords.loc[idx, 'Dim2'], 
                label=serotype, alpha=0.6, color=color_map.get(serotype, '#000000'))
    
# Draw ellipses for each region
regions = famd_row_coords['Region'].unique()
for region in regions:
    region_idx = famd_row_coords['Region'] == region
    confidence_ellipse(famd_row_coords.loc[region_idx, 'Dim1'], 
                       famd_row_coords.loc[region_idx, 'Dim2'], 
                       ax, edgecolor=region_color_map.get(region, '#000000'), alpha=0.3, linewidth=3, facecolor='none')

plt.title('PCA of Primer Combinations by Serotype')
plt.xlabel(f'Component 1 (Variance Explained: {explained_variance_ratio[0]*100:.2f}%)')
plt.ylabel(f'Component 2 (Variance Explained: {explained_variance_ratio[1]*100:.2f}%)')
plt.legend(title='Serotype')
plt.grid(True)

# Save the plot as an image file
output_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\PCA_combined_with_serotype_filtered.png'
plt.savefig(output_file, format='png', dpi=300)

plt.show()

# Print head of the results
print(famd_row_coords.head())
print("Column Contributions:")
print(contrib.head())
print("Row Cos2:")
print(row_cos2.head())

# Save cos2 and contrib to CSV files for further analysis
row_cos2.to_csv(r'C:\Users\Nelso\OneDrive\Documents\Thesis\row_cos2_combined.csv', index=True)
contrib.to_csv(r'C:\Users\Nelso\OneDrive\Documents\Thesis\col_contrib_combined.csv', index=True)
