import pandas as pd
from tqdm.auto import tqdm
import sys

def calculate_average_tm(df_metadata):
    df_metadata['avg_tm'] = (df_metadata['Tm_min'] + df_metadata['Tm_max']) / 2

def prefilter_data(merged_df):
    fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')].copy()
    rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')].copy()
    return fwd_primers, rev_primers

def find_first_round_combinations(fwd_primers, rev_primers, outer_min_amplicon_length, outer_max_amplicon_length):
    combinations = []
    for _, fwd_primer in fwd_primers.iterrows():
        matched_revs = rev_primers[rev_primers['Reference'] == fwd_primer['Reference']]
        amplicon_lengths = abs(matched_revs['Start'] - fwd_primer['Start'])
        valid = matched_revs[(amplicon_lengths >= outer_min_amplicon_length) & (amplicon_lengths <= outer_max_amplicon_length)]
        for _, rev_primer in valid.iterrows():
            combinations.append({
                'Forward_Primer': fwd_primer['Primer'],
                'Reverse_Primer': rev_primer['Primer'],
                'Reference': fwd_primer['Reference'],
                'Amplicon_Length': abs(rev_primer['Start'] - fwd_primer['Start']),
                'Round': 'First'
            })
    return pd.DataFrame(combinations)

def find_second_round_combinations(first_round_combinations, new_primers, inner_min_amplicon_length, inner_max_amplicon_length):
    second_round_combinations = []
    for combo in first_round_combinations.itertuples():
        reused_primer = combo.Forward_Primer if combo.Round == 'First' else combo.Reverse_Primer
        matching_new_primers = new_primers[(new_primers['Reference'] == combo.Reference) & (new_primers['Primer'] != reused_primer)]
        for _, new_primer in matching_new_primers.iterrows():
            amplicon_length = abs(new_primer['Start'] - combo.Amplicon_Length)
            if inner_min_amplicon_length <= amplicon_length <= inner_max_amplicon_length:
                second_round_combinations.append({
                    'Reused_Primer': reused_primer,
                    'New_Primer': new_primer['Primer'],
                    'Reference': combo.Reference,
                    'Amplicon_Length': amplicon_length,
                    'Round': 'Second'
                })
    return pd.DataFrame(second_round_combinations)

def generate_semi_nested_primer_combinations(mapping_positions_path, primer_metadata_path, outer_min_amplicon_length, outer_max_amplicon_length, inner_min_amplicon_length, inner_max_amplicon_length):
    try:
        df_mapping = pd.read_csv(mapping_positions_path, sep='\t')
        df_metadata = pd.read_csv(primer_metadata_path, sep='\t')
        calculate_average_tm(df_metadata)
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        fwd_primers, rev_primers = prefilter_data(merged_df)

        first_round_combinations = find_first_round_combinations(fwd_primers, rev_primers, outer_min_amplicon_length, outer_max_amplicon_length)
        second_round_combinations = find_second_round_combinations(first_round_combinations, fwd_primers.append(rev_primers, ignore_index=True), inner_min_amplicon_length, inner_max_amplicon_length)

        all_combinations = pd.concat([first_round_combinations, second_round_combinations]).drop_duplicates()
        
        # Split into forward and reverse reused
        forward_reused = all_combinations[all_combinations['Reused_Primer'].str.endswith('_F')]
        reverse_reused = all_combinations[all_combinations['Reused_Primer'].str.endswith('_R')]

        if not forward_reused.empty:
            forward_reused.to_csv('semi_nested_forward_reused_u.tsv', sep='\t', index=False)
            print("Forward reused primer combinations saved to 'semi_nested_forward_reused.tsv'")
        
        if not reverse_reused.empty:
            reverse_reused.to_csv('semi_nested_reverse_reused_u.tsv', sep='\t', index=False)
            print("Reverse reused primer combinations saved to 'semi_nested_reverse_reused.tsv'")

        if forward_reused.empty and reverse_reused.empty:
            print("No valid semi-nested primer combinations found.")

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python script.py mapping_positions.tsv primer_metadata.tsv outer_min_amplicon_length outer_max_amplicon_length inner_min_amplicon_length inner_max_amplicon_length")
        sys.exit(1)
    else:
        mapping_positions_path = sys.argv[1]
        primer_metadata_path = sys.argv[2]
        outer_min_amplicon_length = int(sys.argv[3])
        outer_max_amplicon_length = int(sys.argv[4])
        inner_min_amplicon_length = int(sys.argv[5])
        inner_max_amplicon_length = int(sys.argv[6])
        generate_semi_nested_primer_combinations(mapping_positions_path, primer_metadata_path, outer_min_amplicon_length, outer_max_amplicon_length, inner_min_amplicon_length, inner_max_amplicon_length)
