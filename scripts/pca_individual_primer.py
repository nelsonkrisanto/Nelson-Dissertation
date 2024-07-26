import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

# Function to draw confidence ellipse
def confidence_ellipse(x, y, ax, n_std=1.96, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.
    """
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

    transf = (transforms.Affine2D()
              .rotate_deg(45)
              .scale(scale_x, scale_y)
              .translate(mean_x, mean_y))

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

# Load the input data
input_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\individual_conv_primer_details.tsv'  # Update with your actual file path
df = pd.read_csv(input_file, delimiter='\t')

# Display the dataframe to check the loaded data
print(df.head())
print("Columns in DataFrame:", df.columns)

# Store the serotype and primer name information for later use
serotype_info = df[['Primer_Name', 'Genotype', 'Region']].copy()
serotype_info.rename(columns={'Genotype': 'Serotype'}, inplace=True)

# Concatenate columns to create unique row keys
df['Row_Key'] = df.apply(lambda row: f"{row['Primer_Name']}_{row['Orientation']}_{row['Region']}", axis=1)
serotype_info['Row_Key'] = df['Row_Key']  # Add Row_Key to serotype_info for later merging
df.set_index('Row_Key', inplace=True)

# Drop the 'Serotype', 'Primer_Name', and 'Sequence' columns and any other leading information columns not needed
df.drop(columns=['Primer_Name', 'Genotype', 'Sequence'], inplace=True)

# Handle missing values by filling them only for numeric columns
numeric_columns = df.select_dtypes(include=['number']).columns
df[numeric_columns] = df[numeric_columns].fillna(df[numeric_columns].mean())

# Convert categorical columns to category dtype and encode them
categorical_columns = df.select_dtypes(include=['object']).columns
df = pd.get_dummies(df, columns=categorical_columns)

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

# Exclude region-related columns and orientation columns from the contribution calculation
excluded_columns = [col for col in df.columns if col.startswith('Region_') or col.startswith('Orientation_Forward') or col.startswith('Orientation_Reverse')]
non_excluded_columns = [col for col in df.columns if col not in excluded_columns]
contrib = pd.DataFrame(contrib.T, index=['contrib_dim1', 'contrib_dim2'], columns=df.columns).T
contrib = contrib.loc[non_excluded_columns]

# Create a DataFrame for the principal components
famd_row_coords = pd.DataFrame(principal_components, columns=['Dim1', 'Dim2'], index=df.index)

# Merge the coordinates with the serotype and primer name information
famd_row_coords = famd_row_coords.merge(serotype_info, left_index=True, right_on='Row_Key')

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

# Use the 'Serotype' column for classification
serotypes = famd_row_coords['Serotype']
unique_serotypes = serotypes.unique()

# Scatter plot with colors
ax = plt.gca()
for serotype in unique_serotypes:
    idx = serotypes == serotype
    plt.scatter(famd_row_coords.loc[idx, 'Dim1'], famd_row_coords.loc[idx, 'Dim2'], 
                label=serotype, alpha=0.6, color=color_map.get(serotype, '#000000'))
    
    # Draw ellipses for each region
    regions = famd_row_coords.loc[idx, 'Region'].unique()
    for region in regions:
        region_idx = famd_row_coords['Region'] == region
        confidence_ellipse(famd_row_coords.loc[region_idx, 'Dim1'], 
                           famd_row_coords.loc[region_idx, 'Dim2'], 
                           ax, edgecolor=region_color_map.get(region, '#000000'), alpha=0.5, linewidth=2, facecolor='none')

plt.title('PCA of Individual Primers by Serotype')
plt.xlabel(f'Component 1 (Variance Explained: {explained_variance_ratio[0]*100:.2f}%)')
plt.ylabel(f'Component 2 (Variance Explained: {explained_variance_ratio[1]*100:.2f}%)')
plt.legend(title='Serotype')
plt.grid(True)

print(famd_row_coords.head())

# Save the plot as an image file
output_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\PCA_FINAL.png'
plt.savefig(output_file, format='png', dpi=300)
plt.show()

# Print cos2 and contrib
print("Row Cos2 (quality of representation):\n", row_cos2.head())
print("Column Contrib (contribution to PCs):\n", contrib.head())

# Save cos2 and contrib to CSV files for further analysis
row_cos2.to_csv(r'C:\Users\Nelso\OneDrive\Documents\Thesis\row_cos2.csv', index=True)
contrib.to_csv(r'C:\Users\Nelso\OneDrive\Documents\Thesis\col_contrib.csv', index=True)
