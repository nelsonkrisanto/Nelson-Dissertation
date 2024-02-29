import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

# Load the Excel file
file_path = '/home/people/23203786/scratch/Nelson-Dissertation/data/primers/'  # Update this with the path to your Excel file
file = 'test_data.csv'
df = pd.read_excel(file_path+file, header=1)  # Assuming the actual data starts from the second row

# Function to calculate melting temperature
def calculate_tm(sequence):
    if pd.isnull(sequence):
        return None
    seq = Seq(sequence)
    return round(mt.Tm_NN(seq), 2)

# Apply the function to calculate Tm for each sequence
df['Tm°'] = df['Sequence'].apply(calculate_tm)

# Save the updated DataFrame back to an Excel file
output_file = 'up_' + file
output_file_path = file_path + 'output'  # You can change this to your preferred output file path
df.to_excel(output_file_path+output_file, index=False)

print(f"Updated Excel file saved as '{output_file_path}'")
