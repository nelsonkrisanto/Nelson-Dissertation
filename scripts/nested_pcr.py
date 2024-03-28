import pandas as pd
import sys

def calculate_average_tm(min_tm, max_tm):
    """Calculate the average melting temperature (Tm) of a primer."""
    return (min_tm + max_tm) / 2

def generate_nested_primer_combinations(mapping_file, metadata_file):
    try:
        # Load the mapping data
        print(f"Loading mapping data from {mapping_file}")
        df_mapping = pd.read_csv(mapping_file, delimiter='\t')

        # Load the primer metadata
        print(f"Loading primer metadata from {metadata_file}")
        df_metadata = pd.read_csv(metadata_file, delimiter='\t')

        # Merge mapping data with primer metadata
        print("Merging dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        print("Merged dataframes successfully.")

        # List to store first-round primer combinations
        primer_combinations_first_round = []

        # Search for valid first-round primer combinations
        print("Searching for first round valid primer combinations")
        for _, row1 in merged_df.iterrows():
            if row1['Primer'].endswith('_F'):  # Select forward primers
                for _, row2 in merged_df.iterrows():
                    if row2['Primer'].endswith('_R') and row1['Reference'] == row2['Reference']:  # Pair with reverse primers from the same reference
                        # Calculate the amplicon length
                        amplicon_length = abs(row2['Start'] - row1['Start'])
                        if amplicon_length > 300:  # Filter based on amplicon length
                            combination_name = f"{row1['Primer']}_{row2['Primer']}"
                            
                            # Calculate average Tm for both primers
                            avg_tm_fwd = calculate_average_tm(row1['Tm_min'], row1['Tm_max'])
                            avg_tm_rev = calculate_average_tm(row2['Tm_min'], row2['Tm_max'])
                            
                            # Check if the Tm difference is within the acceptable range
                            if abs(avg_tm_fwd - avg_tm_rev) <= 5:
                                # Add the combination to the list
                                primer_combinations_first_round.append({
                                    'Forward_Primer': row1['Primer'],
                                    'Reverse_Primer': row2['Primer'],
                                    'Combination_Name': combination_name,
                                    'Reference': row1['Reference'],
                                    'Amplicon_Length_First_Round': amplicon_length
                                })
                                print(f"Found a valid first-round primer combination: {combination_name}")

        # List to store second-round primer combinations
        primer_combinations_second_round = []

        # Search for valid second-round primer combinations within first-round amplicons
        print("Searching for second round valid primer combinations within first-round amplicons")
        for combo in primer_combinations_first_round:
            # Find the start positions of the first-round primers in the merged_df
            first_round_fwd_start = merged_df.loc[merged_df['Primer'] == combo['Forward_Primer'], 'Start'].values[0]
            first_round_rev_start = merged_df.loc[merged_df['Primer'] == combo['Reverse_Primer'], 'Start'].values[0]

            for _, fwd in merged_df.iterrows():
                if fwd['Primer'].endswith('_F') and fwd['Reference'] == combo['Reference']:
                    for _, rev in merged_df.iterrows():
                        if rev['Primer'].endswith('_R') and rev['Reference'] == combo['Reference']:
                            # Ensure the second-round primers are within the first-round amplicon
                            if first_round_fwd_start < fwd['Start'] < rev['Start'] < first_round_rev_start:
                                amplicon_length_second = abs(rev['Start'] - fwd['Start'])
                                if amplicon_length_second < 100:  # Filter based on nested amplicon length
                                    second_round_name = f"{fwd['Primer']}-{rev['Primer']} within {combo['Combination_Name']}"
                                    primer_combinations_second_round.append({
                                        'First_Round_Combination': combo['Combination_Name'],
                                        'Second_Round_Forward_Primer': fwd['Primer'],
                                        'Second_Round_Reverse_Primer': rev['Primer'],
                                        'Second_Round_Combination_Name': second_round_name,
                                        'Amplicon_Length_Second_Round': amplicon_length_second
                                    })
                                    print(f"Found a valid second-round primer combination: {second_round_name}")


        # Save the second-round primer combinations if any were found
        if primer_combinations_second_round:
            try:
                primer_combinations_df = pd.DataFrame(primer_combinations_second_round)
                primer_combinations_df.drop_duplicates(subset='Second_Round_Combination_Name', inplace=True)

                output_file = 'nested_primer_combinations_updated.tsv'
                primer_combinations_df.to_csv(output_file, sep='\t', index=False)
                print(f"Nested primer combinations saved to {output_file}")

            except Exception as e:
                print(f"An error occurred while saving nested primer combinations: {e}")
        else:
            print("No valid nested primer combinations found.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:  # Expecting two arguments now: the script name, mapping file, and metadata file
        print("Usage: python nested_pcr.py mapping_positions.tsv primer_metadata.tsv")
    else:
        mapping_file = sys.argv[1]
        metadata_file = sys.argv[2]
        generate_nested_primer_combinations(mapping_file, metadata_file)