import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from itertools import product

# Load the Excel file
file_path = 'C:/Users/Nelso/OneDrive/Documents/Thesis/data/'
file = 'cleanned_nodupes_v1.xlsx'
df = pd.read_excel(file_path + file, header=0)

# Ambiguous base mappings
ambiguous_bases = {
    'N': ['A', 'T', 'G', 'C'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G']
}

def generate_sequences(sequence):
    sequence = sequence.upper()
    combinations = [ambiguous_bases.get(nuc, [nuc]) for nuc in sequence]
    return [''.join(seq) for seq in product(*combinations)]

def calculate_tm_avg(sequence, tm_type):
    sequences = generate_sequences(sequence)
    tm_values = []
    for seq in sequences:
        try:
            if tm_type == 'min':
                tm_values.append(mt.Tm_Wallace(Seq(seq)))
            elif tm_type == 'max':
                tm_values.append(mt.Tm_NN(Seq(seq)))
        except ValueError as e:
            # Skip sequences that cause errors
            print(f"Skipping sequence {seq} due to error: {e}")
            continue

    if tm_values:  # Check if the list is not empty
        return round(sum(tm_values) / len(tm_values), 2)
    return None


# Apply the function to calculate the average Tm min and Tm max for each sequence
df['Tm min avg (°C)'] = df['Sequence'].apply(lambda x: calculate_tm_avg(x, 'min'))
df['Tm max avg (°C)'] = df['Sequence'].apply(lambda x: calculate_tm_avg(x, 'max'))

# Save the updated DataFrame back to an Excel file
output_file = 'up_' + file
output_file_path = file_path
df.to_excel(output_file_path + output_file, index=False)

print(f"Updated Excel file saved as '{output_file_path + output_file}'")
