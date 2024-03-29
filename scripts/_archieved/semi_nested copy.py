import pandas as pd
from tqdm import tqdm
import sys

def prefilter_data(merged_df):
    """
    Separates the merged DataFrame into forward and reverse primers for easier processing.

    Parameters:
    - merged_df: The DataFrame containing merged mapping and metadata information.

    Returns:
    - fwd_primers: DataFrame containing forward primers.
    - rev_primers: DataFrame containing reverse primers.
    """
    fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')].copy()
    rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')].copy()
    return fwd_primers, rev_primers

def find_semi_nested_combinations(fwd_primers, rev_primers, is_forward_reused, min_length, max_length):
    """
    Identifies valid semi-nested PCR primer combinations based on the reuse of a primer from the first round.

    Parameters:
    - fwd_primers: DataFrame containing forward primers.
    - rev_primers: DataFrame containing reverse primers.
    - is_forward_reused: Boolean indicating if a forward primer is reused in the second round.
    - min_length: Minimum length for a valid semi-nested PCR amplicon.
    - max_length: Maximum length for a valid semi-nested PCR amplicon.

    Returns:
    - combinations: A list of dictionaries, each representing a valid primer combination.
    """
    combinations = []
    reused_primers, new_primers = (fwd_primers, rev_primers) if is_forward_reused else (rev_primers, fwd_primers)

    for _, reused_primer in tqdm(reused_primers.iterrows(), total=reused_primers.shape[0], desc="Processing Reused Primers"):
        matching_new_primers = new_primers[new_primers['Reference'] == reused_primer['Reference']]
        
        for _, new_primer in matching_new_primers.iterrows():
            amplicon_length = abs(new_primer['Start'] - reused_primer['Start'])
            if min_length <= amplicon_length <= max_length:
                combinations.append({
                    'Reused_Primer': reused_primer['Primer'],
                    'New_Primer': new_primer['Primer'],
                    'Combination_Name': f"{reused_primer['Primer']}_{new_primer['Primer']}",
                    'Reference': reused_primer['Reference'],
                    'Amplicon_Length': amplicon_length
                })
    
    return combinations

def generate_semi_nested_primer_combinations(mapping_positions_path, primer_metadata_path, is_forward_reused):
    """
    Main function to generate semi-nested PCR primer combinations.

    Parameters:
    - mapping_positions_path: Path to the TSV file containing mapping positions.
    - primer_metadata_path: Path to the TSV file containing primer metadata.
    - is_forward_reused: Boolean indicating if a forward primer is reused in the semi-nested PCR.
    """
    try:
        df_mapping = pd.read_csv(mapping_positions_path, sep='\t')
        df_metadata = pd.read_csv(primer_metadata_path, sep='\t')

        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        fwd_primers, rev_primers = prefilter_data(merged_df)

        semi_nested_combinations = find_semi_nested_combinations(fwd_primers, rev_primers, is_forward_reused, 100, 500)
        print(f"Found {len(semi_nested_combinations)} valid semi-nested combinations.")

        if semi_nested_combinations:
            result_df = pd.DataFrame(semi_nested_combinations)
            result_df.drop_duplicates(subset='Combination_Name', inplace=True)
            output_file = 'semi_nested_combinations.tsv'
            result_df.to_csv(output_file, sep='\t', index=False)
            print(f"Semi-nested primer combinations saved to {output_file}")
        else:
            print("No valid semi-nested primer combinations found.")

    except FileNotFoundError as e:
        print(f"File not found: {e}")
        sys.exit(1)
    except pd.errors.EmptyDataError as e:
        print(f"No data: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py mapping_positions.tsv primer_metadata.tsv forward_or_reverse")
        sys.exit(1)
    else:
        mapping_positions_path = sys.argv[1]
        primer_metadata_path = sys.argv[2]
        # Validate the third argument to ensure it's either 'forward' or 'reverse'
        if sys.argv[3].lower() not in ['forward', 'reverse']:
            print("The third argument must be 'forward' or 'reverse'.")
            sys.exit(1)
        is_forward_reused = sys.argv[3].lower() == 'forward'
        generate_semi_nested_primer_combinations(mapping_positions_path, primer_metadata_path, is_forward_reused)
