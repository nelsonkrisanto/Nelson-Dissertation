import pandas as pd
from tqdm import tqdm
import sys

def calculate_average_tm(df_metadata):
    """Calculate and add average Tm directly to the metadata DataFrame."""
    df_metadata['avg_tm'] = (df_metadata['Tm_min'] + df_metadata['Tm_max']) / 2

def prefilter_data(merged_df):
    """Pre-filter the DataFrame to separate forward and reverse primers."""
    fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')].copy()
    rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')].copy()
    return fwd_primers, rev_primers

def find_combinations(fwd_primers, rev_primers, min_length, max_length):
    """Find valid primer combinations using vectorized operations."""
    combinations = []
    for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc="Processing Primers"):
        relevant_revs = rev_primers[(rev_primers['Reference'] == fwd['Reference']) & (rev_primers['Start'] > fwd['Start'])]
        relevant_revs['Amplicon_Length'] = relevant_revs['Start'] - fwd['Start']
        
        valid_revs = relevant_revs[(relevant_revs['Amplicon_Length'] >= min_length) & 
                                   (relevant_revs['Amplicon_Length'] <= max_length) & 
                                   (abs(relevant_revs['avg_tm'] - fwd['avg_tm']) <= 5)]
        
        for _, rev in valid_revs.iterrows():
            combinations.append({
                'Forward_Primer': fwd['Primer'],
                'Reverse_Primer': rev['Primer'],
                'Combination_Name': f"{fwd['Primer']}_{rev['Primer']}",
                'Reference': fwd['Reference'],
                'Amplicon_Length': rev['Amplicon_Length']
            })

    return combinations

def generate_nested_primer_combinations(mapping_file, metadata_file):
    try:
        df_mapping = pd.read_csv(mapping_file, delimiter='\t')
        df_metadata = pd.read_csv(metadata_file, delimiter='\t')
        calculate_average_tm(df_metadata)

        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        fwd_primers, rev_primers = prefilter_data(merged_df)

        first_round_combinations = find_combinations(fwd_primers, rev_primers, 300, float('inf'))
        print(f"Found {len(first_round_combinations)} valid first-round primer combinations.")

        second_round_combinations = find_combinations(fwd_primers, rev_primers, 0, 100)
        print(f"Found {len(second_round_combinations)} valid second-round primer combinations.")

        if second_round_combinations:
            primer_combinations_df = pd.DataFrame(second_round_combinations)
            if 'Combination_Name' in primer_combinations_df.columns:
                primer_combinations_df.drop_duplicates(subset='Combination_Name', inplace=True)
            else:
                print("Column 'Combination_Name' does not exist in DataFrame.")

            output_file = 'nested_primer_combinations_optimized.tsv'
            primer_combinations_df.to_csv(output_file, sep='\t', index=False)
            print(f"Nested primer combinations saved to {output_file}")
        else:
            print("No valid nested primer combinations found.")

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py mapping_positions.tsv primer_metadata.tsv")
        sys.exit(1)
    else:
        mapping_file, metadata_file = sys.argv[1:3]
        generate_nested_primer_combinations(mapping_file, metadata_file)
