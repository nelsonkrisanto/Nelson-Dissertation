import pandas as pd
from tqdm import tqdm
import sys

def calculate_average_tm(df_metadata):
    """
    Calculate and add the average melting temperature (Tm) to the primer metadata DataFrame.
    """
    df_metadata['avg_tm'] = (df_metadata['Tm_min'] + df_metadata['Tm_max']) / 2

def prefilter_data(merged_df):
    """
    Separate the DataFrame into forward and reverse primers for easier processing.
    """
    fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')].copy()
    rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')].copy()
    return fwd_primers, rev_primers

def find_combinations(fwd_primers, rev_primers, min_length, max_length, min_primer_distance):
    """
    Identify valid primer combinations considering amplicon length and minimum primer distance.
    """
    combinations = []
    for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc="Processing Forward Primers"):
        matching_revs = rev_primers[(rev_primers['Reference'] == fwd['Reference']) &
                                    (rev_primers['Start'] > fwd['Start'] + min_primer_distance)]
        matching_revs['Amplicon_Length'] = matching_revs['Start'] - fwd['Start']
        
        valid_revs = matching_revs[(matching_revs['Amplicon_Length'] >= min_length) &
                                   (matching_revs['Amplicon_Length'] <= max_length) &
                                   (abs(matching_revs['avg_tm'] - fwd['avg_tm']) <= 10)]
        
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
    """
    Main function to generate nested primer combinations from mapping and metadata files.
    """
    try:
        df_mapping = pd.read_csv(mapping_file, delimiter='\t')
        df_metadata = pd.read_csv(metadata_file, delimiter='\t')
        calculate_average_tm(df_metadata)

        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        fwd_primers, rev_primers = prefilter_data(merged_df)

        # Define the minimum distance between primer start positions and the amplicon length criteria
        min_primer_distance = 100  # Minimum distance between primer start positions
        min_amplicon_length = 100  # Minimum amplicon length for a valid combination
        max_amplicon_length = 500  # Maximum amplicon length for a valid combination

        # Find valid primer combinations for the first and second rounds with specified criteria
        first_round_combinations = find_combinations(fwd_primers, rev_primers, min_amplicon_length, max_amplicon_length, min_primer_distance)
        second_round_combinations = find_combinations(fwd_primers, rev_primers, min_amplicon_length, max_amplicon_length, min_primer_distance)

        # Save the valid primer combinations to a TSV file
        if first_round_combinations or second_round_combinations:
            all_combinations = first_round_combinations + second_round_combinations
            primer_combinations_df = pd.DataFrame(all_combinations)
            primer_combinations_df.drop_duplicates(subset=['Combination_Name', 'Amplicon_Length', 'Reference'], inplace=True)
            output_file = 'nested_primer_combinations.tsv'
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
