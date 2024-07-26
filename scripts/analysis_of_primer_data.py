import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS

# Load the data from the Excel file
file_path = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\all_cleanned_v1_update.xlsx'
df = pd.read_excel(file_path, sheet_name='Sheet1')  # Update with your actual file path

# Convert 'Year_cleanned' to datetime and extract the year
df['Year_cleanned'] = pd.to_datetime(df['Year_cleanned'], errors='coerce').dt.year

# Drop rows where 'Year_cleanned' is NaN after conversion
df = df.dropna(subset=['Year_cleanned'])

# Ensure 'Year_cleanned' is an integer type
df['Year_cleanned'] = df['Year_cleanned'].astype(int)

# Filter out any year before 1995
df = df[df['Year_cleanned'] >= 1995]

# Standardize the 'Genotype x' column to 'Serotype'
df['Genotype x'] = df['Genotype x'].str.replace('-', '').str.strip()
df.rename(columns={'Genotype x': 'Serotype'}, inplace=True)

# Explode the 'Serotype' column to handle multiple serotypes per row
df['Serotype'] = df['Serotype'].str.split(', ')
df = df.explode('Serotype')

# Group by year and serotype to get counts
year_serotype_counts = df.groupby(['Year_cleanned', 'Serotype']).size().unstack(fill_value=0)

# Save the year-serotype counts to a file
year_serotype_counts.to_csv(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\year_serotype_counts.csv')

# Extract the peak years and counts for each serotype
peak_years = {
    'DENV1': [1999, 2008, 2010, 2018, 2023],
    'DENV2': [2005, 2008, 2010, 2017, 2023],
    'DENV3': [2000, 2003, 2005, 2010, 2018, 2023],
    'DENV4': [2005, 2010, 2014, 2018, 2023],
    'ALL': [2010, 2018, 2023]
}

# Extract the counts for each peak year
peak_counts = {serotype: {year: year_serotype_counts.at[year, serotype] if year in year_serotype_counts.index else 0 for year in years} for serotype, years in peak_years.items()}

# Save the peak counts to a file
with open(r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\peak_counts.txt', 'w') as file:
    for serotype, counts in peak_counts.items():
        file.write(f"{serotype}: {counts}\n")

# Plot settings
sns.set(style="whitegrid")

# Create a bar plot to show the number of primers by year and serotype
plt.figure(figsize=(14, 8))
sns.countplot(data=df, x='Year_cleanned', hue='Serotype', palette='viridis')
plt.title('Number of Primers by Year and Serotype')
plt.xlabel('Year')
plt.ylabel('Count')
plt.legend(title='Serotype')
plt.xticks(rotation=45)
plt.tight_layout()

# Save the plot as an image file
output_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\primers_by_year_and_serotype.png'
plt.savefig(output_file, format='png', dpi=300)
plt.close()

# Create a bar plot to show the number of primers by year and orientation
plt.figure(figsize=(14, 8))
sns.countplot(data=df, x='Year_cleanned', hue='Orientation', palette='viridis')
plt.title('Number of Primers by Year and Orientation')
plt.xlabel('Year')
plt.ylabel('Count')
plt.legend(title='Orientation')
plt.xticks(rotation=45)
plt.tight_layout()

# Save the plot as an image file
output_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\primers_by_year_and_orientation.png'
plt.savefig(output_file, format='png', dpi=300)
plt.close()

# Summary Table
# Total number of primers
total_primers = df['Primer Name'].nunique()

# Total number of publications (titles)
total_publications = df['Title'].nunique()

# Number of duplicates
total_duplicates = df.duplicated().sum()

# Number of primers by orientation
primers_by_orientation = df['Orientation'].value_counts().to_dict()

# Number of primers by serotype
primers_by_serotype = df['Serotype'].value_counts().to_dict()

# Summary table
summary = {
    'Total Primers': total_primers,
    'Total Publications': total_publications,
    'Total Duplicates': total_duplicates,
    'Primers by Orientation': primers_by_orientation,
    'Primers by Serotype': primers_by_serotype
}

# Convert summary to DataFrame for better display
summary_df = pd.DataFrame(list(summary.items()), columns=['Metric', 'Count'])

# Save the summary to a CSV file
summary_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\summary_table.csv'
summary_df.to_csv(summary_file, index=False)

# Additional Summary Information
# Descriptive statistics
descriptive_stats = df.describe(include='all')

# Save descriptive statistics to a CSV file
descriptive_stats_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\descriptive_statistics.csv'
descriptive_stats.to_csv(descriptive_stats_file)

# Primer length statistics
df['Primer Length'] = df['Sequence'].str.len()
primer_length_stats = df['Primer Length'].describe()

# Save primer length statistics to a CSV file
primer_length_stats_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\primer_length_statistics.csv'
primer_length_stats.to_csv(primer_length_stats_file)

# Primer length distribution plot
plt.figure(figsize=(14, 8))
sns.histplot(df['Primer Length'], bins=30, kde=True, color='purple')
plt.title('Distribution of Primer Lengths')
plt.xlabel('Primer Length')
plt.ylabel('Frequency')
plt.tight_layout()

# Save the plot as an image file
output_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\primer_length_distribution.png'
plt.savefig(output_file, format='png', dpi=300)
plt.close()

# Statistical Analysis
yearly_data = df.groupby('Year_cleanned').agg(
    Num_Primers=('Primer Name', 'nunique'),
    Num_Publications=('Title', 'nunique')
).reset_index()

# Calculate weighted number of primers
yearly_data['Weighted_Primers'] = yearly_data['Num_Primers'] / yearly_data['Num_Publications']

# Regression for number of primers
X_primers = sm.add_constant(yearly_data['Year_cleanned'])
y_primers = yearly_data['Num_Primers']
model_primers = OLS(y_primers, X_primers).fit()

# Regression for number of publications
X_publications = sm.add_constant(yearly_data['Year_cleanned'])
y_publications = yearly_data['Num_Publications']
model_publications = OLS(y_publications, X_publications).fit()

# Regression for weighted number of primers
X_weighted_primers = sm.add_constant(yearly_data['Year_cleanned'])
y_weighted_primers = yearly_data['Weighted_Primers']
model_weighted_primers = OLS(y_weighted_primers, X_weighted_primers).fit()

# Save statistical results to a text file
stats_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\statistics_results.txt'
with open(stats_file, 'w') as f:
    f.write("OLS Regression Results for Number of Primers Over the Years:\n")
    f.write(model_primers.summary().as_text())
    f.write("\n\n")
    f.write("OLS Regression Results for Number of Publications Over the Years:\n")
    f.write(model_publications.summary().as_text())
    f.write("\n\n")
    f.write("OLS Regression Results for Weighted Number of Primers Over the Years:\n")
    f.write(model_weighted_primers.summary().as_text())
    f.write("\n\n")
    f.write("Descriptive Statistics:\n")
    f.write(descriptive_stats.to_string())
    f.write("\n\n")
    f.write("Primer Length Statistics:\n")
    f.write(primer_length_stats.to_string())

# Plot the number of primers and publications over the years
fig, ax1 = plt.subplots(figsize=(14, 8))

color = 'tab:blue'
ax1.set_xlabel('Year')
ax1.set_ylabel('Number of Primers', color=color)
ax1.plot(yearly_data['Year_cleanned'], yearly_data['Num_Primers'], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel('Number of Publications', color=color)
ax2.plot(yearly_data['Year_cleanned'], yearly_data['Num_Publications'], color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.title('Number of Primers and Publications Over the Years')
fig.tight_layout()

# Save the plot as an image file
output_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\num_primers_and_publications.png'
plt.savefig(output_file, format='png', dpi=300)
plt.close()

# Plot the weighted number of primers over the years
plt.figure(figsize=(14, 8))
sns.lineplot(data=yearly_data, x='Year_cleanned', y='Weighted_Primers', marker='o', color='green')
plt.title('Weighted Number of Primers Over the Years')
plt.xlabel('Year')
plt.ylabel('Weighted Number of Primers')
plt.tight_layout()

# Save the plot as an image file
output_file = r'C:\Users\Nelso\OneDrive\Documents\Thesis\data\weighted_num_primers.png'
plt.savefig(output_file, format='png', dpi=300)
plt.close()

# Correlation Analysis
correlation, p_value = yearly_data[['Num_Primers', 'Num_Publications']].corr().iloc[0, 1], \
                       yearly_data[['Num_Primers', 'Num_Publications']].corr().iloc[0, 1]

# Save correlation results to a text file
with open(stats_file, 'a') as f:
    f.write(f"Correlation between number of primers and number of publications: {correlation}, P-value: {p_value}")
    f.write("\nSummary of Statistical Analysis:\n")
    f.write(f"Slope for number of primers: {model_primers.params['Year_cleanned']}, P-value: {model_primers.pvalues['Year_cleanned']}\n")
    if model_primers.pvalues['Year_cleanned'] < 0.05:
        f.write("There is a statistically significant increase in the number of primers over the years.\n")
    f.write(f"Slope for number of publications: {model_publications.params['Year_cleanned']}, P-value: {model_publications.pvalues['Year_cleanned']}\n")
    if model_publications.pvalues['Year_cleanned'] < 0.05:
        f.write("There is a statistically significant increase in the number of publications over the years.\n")
    f.write(f"Slope for weighted primers: {model_weighted_primers.params['Year_cleanned']}, P-value: {model_weighted_primers.pvalues['Year_cleanned']}\n")
    if model_weighted_primers.pvalues['Year_cleanned'] < 0.05:
        f.write("There is a statistically significant increase in the weighted number of primers over the years.\n")
    f.write(f"Correlation between number of primers and number of publications: {correlation}, P-value: {p_value}\n")

# Print summary of statistical analysis
print("Summary of Statistical Analysis:")
print(f"Slope for number of primers: {model_primers.params['Year_cleanned']}, P-value: {model_primers.pvalues['Year_cleanned']}")
if model_primers.pvalues['Year_cleanned'] < 0.05:
    print("There is a statistically significant increase in the number of primers over the years.")

print(f"Slope for number of publications: {model_publications.params['Year_cleanned']}, P-value: {model_publications.pvalues['Year_cleanned']}")
if model_publications.pvalues['Year_cleanned'] < 0.05:
    print("There is a statistically significant increase in the number of publications over the years.")

print(f"Slope for weighted primers: {model_weighted_primers.params['Year_cleanned']}, P-value: {model_weighted_primers.pvalues['Year_cleanned']}")
if model_weighted_primers.pvalues['Year_cleanned'] < 0.05:
    print("There is a statistically significant increase in the weighted number of primers over the years.")

print(f"Correlation between number of primers and number of publications: {correlation}, P-value: {p_value}")
