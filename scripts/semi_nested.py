import pandas as pd
from tqdm import tqdm
import sys

def calculate_average_tm(df_metadata):
    """
    Calculates the average melting temperature (Tm) for each primer and adds it as a new column in the DataFrame.
    
    Parameters:
    - df_metadata: DataFrame containing the primer metadata, which includes the minimum and maximum Tm values for each primer.
    """
    # The average Tm is calculated as the mean of the Tm_min and Tm_max values for each primer.
    df_metadata['avg_tm'] = (df_metadata['Tm_min'] + df_metadata['Tm_max']) / 2

def prefilter_data(merged_df):
    """
    Separates the DataFrame containing merged mapping and metadata information into two DataFrames for forward and reverse primers.
    
    Parameters:
    - merged_df: DataFrame resulting from the merging of mapping positions and primer metadata.
    
    Returns:
    - fwd_primers: DataFrame containing only forward primers.
    - rev_primers: DataFrame containing only reverse primers.
    """
    # Filter out forward primers into a separate DataFrame.
    fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')].copy()
    # Filter out reverse primers into another DataFrame.
    rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')].copy()
    return fwd_primers, rev_primers

def find_semi_nested_combinations(fwd_primers, rev_primers, is_forward_reused, min_amplicon_length, max_amplicon_length, min_primer_distance):
    """
    Identifies valid primer combinations for semi-nested PCR, taking into account the amplicon length, 
    the minimum distance between primer start positions, and the compatibility of average melting temperatures.
    
    Parameters:
    - fwd_primers: DataFrame of forward primers.
    - rev_primers: DataFrame of reverse primers.
    - is_forward_reused: Boolean indicating whether the forward primer is reused in the semi-nested PCR.
    - min_amplicon_length: Minimum length required for a valid amplicon.
    - max_amplicon_length: Maximum length allowed for a valid amplicon.
    - min_primer_distance: Minimum distance required between the start positions of the forward and reverse primers.
    
    Returns:
    - combinations: A list of dictionaries, each representing a valid primer combination with details.
    """
    combinations = []
    # Choose which set of primers will be reused based on is_forward_reused flag.
    reused_primers, new_primers = (fwd_primers, rev_primers) if is_forward_reused else (rev_primers, fwd_primers)

    # Iterate through reused primers to find compatible new primers.
    for _, reused_primer in tqdm(reused_primers.iterrows(), total=reused_primers.shape[0], desc="Processing Reused Primers"):
        # Filter new primers by matching reference and ensuring the required distance from the reused primer.
        matching_new_primers = new_primers[(new_primers['Reference'] == reused_primer['Reference'])]

        for _, new_primer in matching_new_primers.iterrows():
            # Calculate the amplicon length and check if it falls within the specified range.
            amplicon_length = abs(new_primer['Start'] - reused_primer['Start'])
            if (min_amplicon_length <= amplicon_length <= max_amplicon_length) and \
               (abs(new_primer['avg_tm'] - reused_primer['avg_tm']) <= 10):  # Ensure Tm compatibility within 10 units.
                combinations.append({
                    'Reused_Primer': reused_primer['Primer'],
                    'New_Primer': new_primer['Primer'],
                    'Combination_Name': f"{reused_primer['Primer']}_{new_primer['Primer']}",
                    'Reference': reused_primer['Reference'],
                    'Amplicon_Length': amplicon_length
                })

    return combinations

def generate_semi_nested_primer_combinations(mapping_positions_path, primer_metadata_path):
    """
    Main function to orchestrate the generation of semi-nested PCR primer combinations from mapping positions and metadata files.
    
    Parameters:
    - mapping_positions_path: File path to the TSV file containing mapping positions.
    - primer_metadata_path: File path to the TSV file containing primer metadata.
    """
    try:
        # Load the mapping positions and primer metadata into DataFrames.
        df_mapping = pd.read_csv(mapping_positions_path, sep='\t')
        df_metadata = pd.read_csv(primer_metadata_path, sep='\t')

        # Calculate the average Tm for each primer.
        calculate_average_tm(df_metadata)

        # Merge the mapping data with the primer metadata.
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')

        # Separate the merged DataFrame into forward and reverse primers.
        fwd_primers, rev_primers = prefilter_data(merged_df)

        # Define criteria for amplicon length and minimum primer distance.
        min_amplicon_length = 100  # Minimum amplicon length for valid combinations.
        max_amplicon_length = 500  # Maximum amplicon length for valid combinations.
        min_primer_distance = 100  # Minimum distance between primer start positions.

        # Find valid primer combinations for both forward and reverse reused primers.
        forward_combinations = find_semi_nested_combinations(fwd_primers, rev_primers, True, min_amplicon_length, max_amplicon_length, min_primer_distance)
        reverse_combinations = find_semi_nested_combinations(fwd_primers, rev_primers, False, min_amplicon_length, max_amplicon_length, min_primer_distance)

        # Save the combinations to separate files based on whether a forward or reverse primer is reused.
        if forward_combinations:
            pd.DataFrame(forward_combinations).drop_duplicates().to_csv('semi_nested_combinations_forward.tsv', sep='\t', index=False)
            print("Forward primer combinations saved to 'semi_nested_combinations_forward.tsv'")

        if reverse_combinations:
            pd.DataFrame(reverse_combinations).drop_duplicates().to_csv('semi_nested_combinations_reverse.tsv', sep='\t', index=False)
            print("Reverse primer combinations saved to 'semi_nested_combinations_reverse.tsv'")

        if not forward_combinations and not reverse_combinations:
            print("No valid semi-nested primer combinations found.")

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py mapping_positions.tsv primer_metadata.tsv")
        sys.exit(1)
    else:
        mapping_positions_path = sys.argv[1]
        primer_metadata_path = sys.argv[2]
        generate_semi_nested_primer_combinations(mapping_positions_path, primer_metadata_path)
