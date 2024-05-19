import pandas as pd
from tqdm import tqdm
import sys

def calculate_average_tm(df_metadata):
    df_metadata['avg_tm'] = (df_metadata['Tm_min'] + df_metadata['Tm_max']) / 2

def prefilter_data(merged_df):
    fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')].copy()
    rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')].copy()
    return fwd_primers, rev_primers

def find_combinations(fwd_primers, rev_primers, min_length, max_length, min_primer_distance, round_type='first'):
    # Create a DataFrame to store the combinations
    combinations = pd.DataFrame()

    # Iterate over forward primers
    for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc=f"Processing {round_type.capitalize()} Round Forward Primers"):
        # Filter reverse primers based on reference and position
        matching_revs = rev_primers[(rev_primers['Reference'] == fwd['Reference']) & (rev_primers['Start'] > fwd['Start'] + min_primer_distance)]
        matching_revs['Amplicon_Length'] = matching_revs['Start'] - fwd['Start']

        # Filter reverse primers based on amplicon length and Tm difference
        valid_revs = matching_revs[(matching_revs['Amplicon_Length'] >= min_length) & (matching_revs['Amplicon_Length'] <= max_length) & (abs(matching_revs['avg_tm'] - fwd['avg_tm']) <= 10)]
        
        # Create a temporary DataFrame with the current forward primer and all valid reverse primers
        temp_df = pd.DataFrame({
            'Forward_Primer': fwd['Primer'],
            'Reverse_Primer': valid_revs['Primer'],
            'Combination_Name': fwd['Primer'] + '_' + valid_revs['Primer'],
            'Reference': fwd['Reference'],
            'Amplicon_Length': valid_revs['Amplicon_Length'],
            'Round_Type': round_type
        })

        # Append the temporary DataFrame to the combinations DataFrame
        combinations = pd.concat([combinations, temp_df], ignore_index=True)

    return combinations

def generate_nested_primer_combinations(mapping_file, metadata_file, outer_min_len, outer_max_len, inner_min_len, inner_max_len):
    try:
        df_mapping = pd.read_csv(mapping_file, delimiter='\t')
        df_metadata = pd.read_csv(metadata_file, delimiter='\t')
        calculate_average_tm(df_metadata)

        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        fwd_primers, rev_primers = prefilter_data(merged_df)

        min_primer_distance = 50  # Minimum distance between primer start positions for both rounds

        # Find valid primer combinations for the first round with specified outer amplicon length criteria
        first_round_combinations = find_combinations(fwd_primers, rev_primers, outer_min_len, outer_max_len, min_primer_distance)

        # For the second round, we find valid primer combinations within the first-round amplicons using the inner amplicon length criteria
        second_round_combinations = pd.DataFrame()
        for combo in first_round_combinations.itertuples(index=False):
            specific_fwd = fwd_primers[fwd_primers['Primer'] == combo.Forward_Primer].copy()
            specific_rev = rev_primers[rev_primers['Primer'] == combo.Reverse_Primer].copy()
            specific_fwd['First_Round_End'] = specific_rev['End'].values[0]  # End position of the first-round amplicon
            temp_combos = find_combinations(specific_fwd, specific_rev, inner_min_len, inner_max_len, 0, 'second')
            second_round_combinations = pd.concat([second_round_combinations, temp_combos], ignore_index=True)

        all_combinations = pd.concat([first_round_combinations, second_round_combinations], ignore_index=True)
        if not all_combinations.empty:
            all_combinations.drop_duplicates(subset=['Combination_Name', 'Amplicon_Length', 'Reference'], inplace=True)
            output_file = 'nested_primer_combinations_ops_u.tsv'
            all_combinations.to_csv(output_file, sep='\t', index=False)
            print(f"Nested primer combinations saved to {output_file}")
        else:
            print("No valid nested primer combinations found.")
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python script.py mapping_positions.tsv primer_metadata.tsv outer_min_len outer_max_len inner_min_len inner_max_len")
        sys.exit(1)
    else:
        mapping_file, metadata_file, outer_min_len, outer_max_len, inner_min_len, inner_max_len = sys.argv[1:7]
        generate_nested_primer_combinations(mapping_file, metadata_file, int(outer_min_len), int(outer_max_len), int(inner_min_len), int(inner_max_len))
