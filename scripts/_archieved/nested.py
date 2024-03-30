import pandas as pd
from tqdm import tqdm
import sys

def calculate_average_tm(df_metadata):
    """
    Adds a column for the average melting temperature (Tm) to the primer metadata DataFrame.
    
    Parameters:
    - df_metadata: DataFrame containing primer metadata, including minimum and maximum Tm values.
    """
    # Calculate average Tm and store it in a new column 'avg_tm'
    df_metadata['avg_tm'] = (df_metadata['Tm_min'] + df_metadata['Tm_max']) / 2

def prefilter_data(merged_df):
    """
    Separates the merged DataFrame into forward and reverse primers for easier processing.
    
    Parameters:
    - merged_df: DataFrame resulting from merging mapping and metadata information on primers.
    
    Returns:
    - Two DataFrames: one for forward primers and another for reverse primers.
    """
    # Filter forward primers and create a copy to avoid SettingWithCopyWarning
    fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')].copy()
    # Filter reverse primers and create a copy
    rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')].copy()
    return fwd_primers, rev_primers

def find_combinations(fwd_primers, rev_primers, min_length, max_length):
    """
    Identifies valid primer combinations based on amplicon length and Tm compatibility.
    
    Parameters:
    - fwd_primers: DataFrame of forward primers.
    - rev_primers: DataFrame of reverse primers.
    - min_length: Minimum acceptable amplicon length for a valid combination.
    - max_length: Maximum acceptable amplicon length for a valid combination.
    
    Returns:
    - A list of dictionaries, each representing a valid primer combination.
    """
    combinations = []  # Initialize an empty list to store valid combinations.
    # Iterate through forward primers with progress tracking
    for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc="Processing Primers"):
        # Find reverse primers from the same reference and within the specified amplicon length range
        relevant_revs = rev_primers[(rev_primers['Reference'] == fwd['Reference']) & 
                                    (rev_primers['Start'] > fwd['Start'])]
        relevant_revs['Amplicon_Length'] = relevant_revs['Start'] - fwd['Start']
        
        # Filter reverse primers based on amplicon length and average Tm difference
        valid_revs = relevant_revs[(relevant_revs['Amplicon_Length'] >= min_length) & 
                                   (relevant_revs['Amplicon_Length'] <= max_length) & 
                                   (abs(relevant_revs['avg_tm'] - fwd['avg_tm']) <= 10)]
        
        # For each valid reverse primer, create a dictionary with combination details and add to the list
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
    
    Parameters:
    - mapping_file: Path to the TSV file containing mapping positions.
    - metadata_file: Path to the TSV file containing primer metadata.
    """
    try:
        # Load mapping positions and primer metadata from provided files
        df_mapping = pd.read_csv(mapping_file, delimiter='\t')
        df_metadata = pd.read_csv(metadata_file, delimiter='\t')
        
        # Calculate average Tm for each primer in metadata
        calculate_average_tm(df_metadata)

        # Merge mapping data with primer metadata on 'Primer' column
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        
        # Separate merged data into forward and reverse primers
        fwd_primers, rev_primers = prefilter_data(merged_df)

        # Find valid first-round and second-round primer combinations with specified amplicon length criteria
        first_round_combinations = find_combinations(fwd_primers, rev_primers, 300, float('inf'))
        print(f"Found {len(first_round_combinations)} valid first-round primer combinations.")

        second_round_combinations = find_combinations(fwd_primers, rev_primers, 0, 100)
        print(f"Found {len(second_round_combinations)} valid second-round primer combinations.")

        # If valid second-round combinations are found, save them to a TSV file
        if second_round_combinations:
            primer_combinations_df = pd.DataFrame(second_round_combinations)
            if 'Combination_Name' in primer_combinations_df.columns and 'Amplicon_Length' in primer_combinations_df.columns:
                primer_combinations_df.drop_duplicates(subset=['Combination_Name', 'Amplicon_Length', 'Reference'], inplace=True)
            else:
                print("Columns 'Combination_Name', 'Amplicon_Length', and 'Reference' do not exist in DataFrame.")


            # Define output file path and save the DataFrame as a TSV file
            output_file = 'nested_combinations.tsv'
            primer_combinations_df.to_csv(output_file, sep='\t', index=False)
            print(f"Nested primer combinations saved to {output_file}")
        else:
            print("No valid nested primer combinations found.")

    except Exception as e:
        # Catch and print any errors that occur during processing
        print(f"An error occurred: {e}")
        sys.exit(1)  # Exit with an error status

if __name__ == "__main__":
    # Check for correct usage and number of command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py mapping_positions.tsv primer_metadata.tsv")
        sys.exit(1)  # Exit with an error status if incorrect usage
    else:
        # Extract file paths from command-line arguments
        mapping_file, metadata_file = sys.argv[1:3]
        # Generate nested primer combinations with the provided files
        generate_nested_primer_combinations(mapping_file, metadata_file)
