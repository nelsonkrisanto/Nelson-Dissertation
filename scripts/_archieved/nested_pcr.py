import pandas as pd
import sys
from tqdm import tqdm

def calculate_average_tm(min_tm, max_tm):
    """
    Calculate the average melting temperature (Tm) of a primer.
    
    Parameters:
    - min_tm: The minimum melting temperature of the primer.
    - max_tm: The maximum melting temperature of the primer.
    
    Returns:
    - The average melting temperature.
    """
    return (min_tm + max_tm) / 2

def generate_nested_primer_combinations(mapping_file, metadata_file):
    """
    Generate nested primer combinations from mapping and metadata files.
    
    Parameters:
    - mapping_file: The path to the CSV file containing mapping positions.
    - metadata_file: The path to the CSV file containing primer metadata.
    """
    try:
        # Load the mapping data from the provided CSV file
        print(f"Loading mapping data from {mapping_file}")
        df_mapping = pd.read_csv(mapping_file, delimiter='\t')

        # Load the primer metadata from the provided CSV file
        print(f"Loading primer metadata from {metadata_file}")
        df_metadata = pd.read_csv(metadata_file, delimiter='\t')

        # Merge the mapping and metadata dataframes based on the 'Primer' column
        print("Merging dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        print("Merged dataframes successfully.")

        # Initialize a list to store valid first-round primer combinations
        primer_combinations_first_round = []

        # Iterate through each primer in the merged dataframe to find valid first-round primer combinations
        print("Searching for first round valid primer combinations")
        for _, row1 in tqdm(merged_df.iterrows(), total=merged_df.shape[0], desc="First Round Progress"):
            if row1['Primer'].endswith('_F'):  # Check if the primer is a forward primer
                for _, row2 in merged_df.iterrows():
                    if row2['Primer'].endswith('_R') and row1['Reference'] == row2['Reference']:  # Check for matching reverse primer and reference
                        amplicon_length = abs(row2['Start'] - row1['Start'])  # Calculate the amplicon length
                        if amplicon_length > 300:  # Filter based on amplicon length criteria for the first round
                            combination_name = f"{row1['Primer']}_{row2['Primer']}"
                            
                            # Calculate average Tm for both primers
                            avg_tm_fwd = calculate_average_tm(row1['Tm_min'], row1['Tm_max'])
                            avg_tm_rev = calculate_average_tm(row2['Tm_min'], row2['Tm_max'])
                            
                            # Check if the Tm difference is within the acceptable range
                            if abs(avg_tm_fwd - avg_tm_rev) <= 5:
                                # Add the valid combination to the list
                                primer_combinations_first_round.append({
                                    'Forward_Primer': row1['Primer'],
                                    'Reverse_Primer': row2['Primer'],
                                    'Combination_Name': combination_name,
                                    'Reference': row1['Reference'],
                                    'Amplicon_Length_First_Round': amplicon_length
                                })

        # Initialize a list to store valid second-round (nested) primer combinations
        primer_combinations_second_round = []

        # Iterate through the valid first-round combinations to find valid second-round primer combinations
        print("Searching for second round valid primer combinations within first-round amplicons")
        for combo in tqdm(primer_combinations_first_round, desc="Second Round Progress"):
            # Extract the start positions of the first-round primers
            first_round_fwd_start = merged_df.loc[merged_df['Primer'] == combo['Forward_Primer'], 'Start'].values[0]
            first_round_rev_start = merged_df.loc[merged_df['Primer'] == combo['Reverse_Primer'], 'Start'].values[0]

            # Iterate through each primer again to find valid second-round combinations
            for _, fwd in merged_df.iterrows():
                if fwd['Primer'].endswith('_F') and fwd['Reference'] == combo['Reference']:
                    for _, rev in merged_df.iterrows():
                        if rev['Primer'].endswith('_R') and rev['Reference'] == combo['Reference']:
                            # Ensure the second-round primers are within the first-round amplicon
                            if first_round_fwd_start < fwd['Start'] < rev['Start'] < first_round_rev_start:
                                amplicon_length_second = abs(rev['Start'] - fwd['Start'])  # Calculate the second-round amplicon length
                                if amplicon_length_second < 100:  # Filter based on amplicon length criteria for the second round
                                    second_round_name = f"{fwd['Primer']}-{rev['Primer']} within {combo['Combination_Name']}"
                                    # Add the valid second-round combination to the list
                                    primer_combinations_second_round.append({
                                        'First_Round_Combination': combo['Combination_Name'],
                                        'Second_Round_Forward_Primer': fwd['Primer'],
                                        'Second_Round_Reverse_Primer': rev['Primer'],
                                        'Second_Round_Combination_Name': second_round_name,
                                        'Amplicon_Length_Second_Round': amplicon_length_second
                                    })

        # Save the valid second-round primer combinations to a CSV file if any were found
        if primer_combinations_second_round:
            primer_combinations_df = pd.DataFrame(primer_combinations_second_round)
            primer_combinations_df.drop_duplicates(subset='Second_Round_Combination_Name', inplace=True)  # Remove any duplicate combinations

            output_file = 'nested_primer_combinations_updated.tsv'  # Specify the output file name
            primer_combinations_df.to_csv(output_file, sep='\t', index=False)  # Save the dataframe to a TSV file
            print(f"Nested primer combinations saved to {output_file}")
        else:
            print("No valid nested primer combinations found.")

    except Exception as e:
        print(f"An error occurred: {e}")  # Print any errors that occur during the execution
        sys.exit(1)  # Exit the script with an error status

if __name__ == "__main__":
    # Check if the correct number of command line arguments were provided
    if len(sys.argv) != 3:
        print("Usage: python nested_pcr.py mapping_positions.tsv primer_metadata.tsv")
        sys.exit(1)  # Exit the script with an error status if incorrect number of arguments
    else:
        mapping_file = sys.argv[1]  # Assign the first argument to mapping_file
        metadata_file = sys.argv[2]  # Assign the second argument to metadata_file
        generate_nested_primer_combinations(mapping_file, metadata_file)  # Call the function with the provided file paths
