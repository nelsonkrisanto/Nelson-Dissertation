import pandas as pd
import sys

def generate_primer_combinations(input_file):
    try:
        print(f"Loading mapping data from {input_file}")
        # Load the dataframe from the input TSV file
        df_mapping = pd.read_csv(input_file, delimiter='\t')

        print("Loading primer metadata from primer_metadata.tsv")
        # Load the primer metadata from the primer_metadata.tsv file
        df_metadata = pd.read_csv('primer_metadata.tsv', delimiter='\t')

        print("Merging dataframes based on 'Primer' column")
        # Merge mapping and metadata dataframes on 'Primer' column
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        print("Merged dataframes successfully.")

        # Initialize a list to store primer combinations
        primer_combinations = []

        print("Searching for valid primer combinations")
        # Iterate through the merged dataframe to find valid primer combinations
        for _, row1 in merged_df.iterrows():
            if row1['Primer'].endswith('-F') and row1['orientation'] == 'forward':
                for _, row2 in merged_df.iterrows():
                    if row2['Primer'].endswith('-R') and row1['Reference'] == row2['Reference'] and row2['Genogroup'] == row1['Genogroup'] and row2['orientation'] == 'reverse':
                        diff = abs(row1['Start'] - row2['Start'])
                        if diff > 150:
                            amplicon_length = row2['Start'] - row1['Start']
                            if 100 < amplicon_length < 1000:  # Filter amplicon length
                                combination_name = f"{row1['Primer']}_{row2['Primer']}"
                                
                                # Calculate absolute and rounded differences
                                tm_max_diff = round(abs(row1['Tm_max'] - row2['Tm_max']))
                                tm_min_diff = round(abs(row1['Tm_min'] - row2['Tm_min']))
                                
                                if tm_max_diff <= 7 and tm_min_diff <= 7:
                                    primer_combinations.append({
                                        'Forward_Primer': row1['Primer'],
                                        'Reverse_Primer': row2['Primer'],
                                        'Combination_Name': combination_name,
                                        'Reference': row1['Reference'],
                                        'Amplicon_Length': amplicon_length,
                                        'Genogroup': row1['Genogroup'],
                                        'Forward_Start': row1['Start'],
                                        'Forward_End': row1['End'],
                                        'Reverse_Start': row2['Start']
                                    })
                                    print(f"Found a valid primer combination: {combination_name}")

        if primer_combinations:
            try:
                # Create a new dataframe from the primer_combinations list
                primer_combinations_df = pd.DataFrame(primer_combinations)

                # Remove duplicate combinations based on 'Combination_Name'
                primer_combinations_df.drop_duplicates(subset='Combination_Name', inplace=True)

                print("Saving primer combinations to primer_combinations.tsv")
                # Save the primer combinations dataframe to primer_combinations.tsv
                output_file = 'primer_combinations.tsv'
                primer_combinations_df.to_csv(output_file, sep='\t', index=False)
                print(f"Primer combinations saved to {output_file}")

                # Create a new dataframe for insilico_pcr_primers.tsv
                insilico_pcr_df = pd.DataFrame(primer_combinations_df, columns=['Forward_Primer', 'Reverse_Primer', 'Combination_Name'])
                insilico_pcr_output_file = 'insilico_pcr_primers.tsv'
                insilico_pcr_df.to_csv(insilico_pcr_output_file, sep='\t', index=False, header=False)
                print(f"insilico_pcr_primers saved to {insilico_pcr_output_file}")
            except Exception as e:
                print(f"An error occurred while saving primer combinations: {e}")
        else:
            print("No valid primer combinations found.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:  # Expecting two arguments: script name and input file
        print("Usage: python combine.py mapping_positions.tsv primer_metadata.tsv")
    else:
        input_file = sys.argv[1]
        generate_primer_combinations(input_file)
