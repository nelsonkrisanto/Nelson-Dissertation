import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

# Load the Excel file
file_path = 'C:/Users/Nelso/OneDrive/Documents/Thesis/data/'  # Update this with the path to your Excel file
file = 'cleanned_nodupes_v1.xlsx'
df = pd.read_excel(file_path + file, header=0)  # Use header=0 if column names are in the first row

# Function to replace ambiguous bases for Tm min calculation (choose weaker bonds)
def replace_ambiguous_min(sequence):
    return sequence.upper().replace('N', 'A').replace('R', 'A').replace('Y', 'T')\
                           .replace('S', 'C').replace('W', 'A').replace('K', 'T')\
                           .replace('M', 'A').replace('B', 'T').replace('D', 'A')\
                           .replace('H', 'A').replace('V', 'A')

# Function to replace ambiguous bases for Tm max calculation (choose stronger bonds)
def replace_ambiguous_max(sequence):
    return sequence.upper().replace('N', 'G').replace('R', 'G').replace('Y', 'C')\
                           .replace('S', 'G').replace('W', 'G').replace('K', 'G')\
                           .replace('M', 'C').replace('B', 'G').replace('D', 'G')\
                           .replace('H', 'C').replace('V', 'G')

# Adjusted function to calculate minimum melting temperature
def calculate_tm_min(sequence):
    if pd.isnull(sequence):
        return None
    seq = Seq(replace_ambiguous_min(sequence))
    return round(mt.Tm_Wallace(seq), 2)

# Adjusted function to calculate maximum melting temperature
def calculate_tm_max(sequence):
    if pd.isnull(sequence):
        return None
    seq = Seq(replace_ambiguous_max(sequence))
    try:
        return round(mt.Tm_NN(seq), 2)
    except ValueError as e:
        print(f"Error processing sequence {sequence}: {e}")
        return None

# Apply the functions to calculate Tm min and Tm max for each sequence
df['Tm min (°C)'] = df['Sequence'].apply(calculate_tm_min)
df['Tm max (°C)'] = df['Sequence'].apply(calculate_tm_max)

# Save the updated DataFrame back to an Excel file
output_file = 'up_' + file
output_file_path = file_path
df.to_excel(output_file_path + output_file, index=False)

print(f"Updated Excel file saved as '{output_file_path + output_file}'")
