import pandas as pd
import sys
from tqdm.auto import tqdm

# Function to calculate the average melting temperature (Tm) for each primer
def calculate_average_tm(df_metadata):
    # This is a vectorized operation, applying the calculation to the entire column at once, which is faster than iterating over rows
    df_metadata['avg_tm'] = (df_metadata['Tm_min'] + df_metadata['Tm_max']) / 2

# Function to separate forward and reverse primers based on a suffix in their names
def prefilter_data(merged_df):
    # Use boolean indexing to filter rows, which is more efficient than using a loop
    fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')]
    rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')]
    return fwd_primers, rev_primers

# Function to find valid primer combinations based on several criteria
def find_combinations(fwd_primers, rev_primers, min_length, max_length, min_primer_distance, round_type='first'):
    # Rename columns in reverse primers DataFrame to avoid column name conflicts after merging
    rev_primers = rev_primers.rename(columns={'Start': 'Rev_Start', 'avg_tm': 'Rev_avg_tm'})
    
    # Merge forward and reverse primers on the 'Reference' column to find matching pairs within the same reference sequence
    merged = fwd_primers.merge(rev_primers, on='Reference')
    
    # Filter merged DataFrame for primer pairs that meet the minimum distance requirement
    merged = merged[merged['Rev_Start'] > merged['Start'] + min_primer_distance]
    
    # Calculate the amplicon length for each primer pair
    merged['Amplicon_Length'] = merged['Rev_Start'] - merged['Start']
    
    # Further filter primer pairs based on amplicon length and Tm difference criteria
    valid_revs = merged[(merged['Amplicon_Length'] >= min_length) & 
                        (merged['Amplicon_Length'] <= max_length) & 
                        (abs(merged['Rev_avg_tm'] - merged['avg_tm']) <= 10)]
    
    # Construct columns for the output DataFrame from the filtered primer pairs
    valid_revs['Forward_Primer'] = valid_revs['Primer_x']
    valid_revs['Reverse_Primer'] = valid_revs['Primer_y']
    valid_revs['Combination_Name'] = valid_revs['Forward_Primer'] + "_" + valid_revs['Reverse_Primer']
    valid_revs['Round_Type'] = round_type
    
    # Select relevant columns for the final output
    combinations = valid_revs[['Forward_Primer', 'Reverse_Primer', 'Combination_Name', 'Reference', 'Amplicon_Length', 'Round_Type']]
    
    # Convert DataFrame to a list of dictionaries for easier handling or further processing
    return combinations.to_dict('records')

# Main function to orchestrate the process of generating nested primer combinations
def generate_nested_primer_combinations(mapping_file, metadata_file, outer_min_len, outer_max_len, inner_min_len, inner_max_len):
    try:
        # Load mapping and metadata information from provided files
        df_mapping = pd.read_csv(mapping_file, delimiter='\t')
        df_metadata = pd.read_csv(metadata_file, delimiter='\t')
        
        # Calculate average Tm for each primer
        calculate_average_tm(df_metadata)

        # Merge mapping and metadata DataFrames on the 'Primer' column
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        
        # Prefilter to separate forward and reverse primers
        fwd_primers, rev_primers = prefilter_data(merged_df)

        # Set minimum distance between primer start positions for both rounds
        min_primer_distance = 100

        # Find valid primer combinations for the first round with specified outer amplicon length criteria
        first_round_combinations = find_combinations(fwd_primers, rev_primers, outer_min_len, outer_max_len, min_primer_distance)

        # For the second round, find valid primer combinations within the first-round amplicons using the inner amplicon length criteria
        second_round_combinations = []
        for combo in tqdm(first_round_combinations, desc='Processing Second Round Combinations'):
            specific_fwd = fwd_primers[fwd_primers['Primer'] == combo['Forward_Primer']]
            specific_rev = rev_primers[rev_primers['Primer'] == combo['Reverse_Primer']]
            
            # Update the end position for the first-round amplicons to limit the search space for the second round
            specific_fwd['First_Round_End'] = specific_rev['End'].values[0]
            
            # Find second roundcombinations += find_combinations(specific_fwd, specific_rev, inner_min_len, inner_max_len, 0, 'second')

        # Concatenate first and second round combinations and remove duplicates
        all_combinations = pd.DataFrame(first_round_combinations + second_round_combinations).drop_duplicates(subset=['Combination_Name', 'Amplicon_Length', 'Reference'])

        # Check if any combinations were found and save to file if so
        if not all_combinations.empty:
            output_file = 'nested_primer_combinations_ops.tsv'
            all_combinations.to_csv(output_file, sep='\t', index=False)
            print(f"Nested primer combinations saved to {output_file}")
        else:
            print("No valid nested primer combinations found.")
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

# Entry point of the script
if __name__ == "__main__":
    # Ensure the correct number of arguments are provided
    if len(sys.argv) != 7:
        print("Usage: python script.py mapping_positions.tsv primer_metadata.tsv outer_min_len outer_max_len inner_min_len inner_max_len")
        sys.exit(1)
    else:
        # Parse arguments and call the main function with the parsed parameters
        mapping_file, metadata_file, outer_min_len, outer_max_len, inner_min_len, inner_max_len = sys.argv[1:7]
        generate_nested_primer_combinations(mapping_file, metadata_file, int(outer_min_len), int(outer_max_len), int(inner_min_len), int(inner_max_len))
