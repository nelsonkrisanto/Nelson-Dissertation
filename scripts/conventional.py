import pandas as pd
import sys
from tqdm import tqdm

def generate_primer_combinations(input_file):
    """
    Generate valid primer combinations based on criteria from input files.

    Parameters:
    - input_file: The file path to the TSV file containing mapping positions.
    """
    try:
        # Load the mapping data from the provided input TSV file.
        print(f"Loading mapping data from {input_file}")
        df_mapping = pd.read_csv(input_file, delimiter='\t')

        # Load primer metadata from a fixed file named 'primer_metadata.tsv'.
        print("Loading primer metadata from primer_metadata.tsv")
        df_metadata = pd.read_csv('primer_metadata.tsv', delimiter='\t')

        # Merge the mapping and metadata dataframes based on the 'Primer' column.
        # This combines related information for each primer from both files into a single DataFrame.
        print("Merging dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        print("Merged dataframes successfully.")

        # Separate the merged DataFrame into two: one for forward primers and one for reverse primers.
        fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')]
        rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')]

        # Initialize a list to store information about valid primer combinations.
        primer_combinations = []

        # Iterate through each forward primer to find matching reverse primers based on specific criteria.
        print("Searching for valid primer combinations")
        for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc="Processing Primers"):
            # Filter reverse primers that are from the same reference as the forward primer
            # and where the start position is appropriately distanced to ensure a valid amplicon length.
            matching_revs = rev_primers[(rev_primers['Reference'] == fwd['Reference']) & 
                                        (rev_primers['Start'] > fwd['Start'] + 150) & 
                                        (rev_primers['Start'] < fwd['Start'] + 1000)]
            
            # For each matching reverse primer, validate the primer combination based on Tm values and amplicon length.
            for _, rev in matching_revs.iterrows():
                combination_name = f"{fwd['Primer']}_{rev['Primer']}"
                amplicon_length = rev['Start'] - fwd['Start']  # Calculate the length of the potential amplicon.
                tm_max_diff = round(abs(fwd['Tm_max'] - rev['Tm_max']))  # Calculate the difference in max Tm values.
                tm_min_diff = round(abs(fwd['Tm_min'] - rev['Tm_min']))  # Calculate the difference in min Tm values.

                # Check if the Tm differences are within the acceptable range.
                if tm_max_diff <= 7 and tm_min_diff <= 7:
                    # If all criteria are met, add the primer combination to the list.
                    primer_combinations.append({
                        'Forward_Primer': fwd['Primer'],
                        'Reverse_Primer': rev['Primer'],
                        'Combination_Name': combination_name,
                        'Reference': fwd['Reference'],
                        'Amplicon_Length': amplicon_length
                    })

        # If any valid primer combinations were found, save them to a TSV file.
        if primer_combinations:
            # Convert the list of primer combinations to a DataFrame.
            primer_combinations_df = pd.DataFrame(primer_combinations)
            # Remove any duplicate combinations based on 'Combination_Name'.
            primer_combinations_df.drop_duplicates(subset='Combination_Name', inplace=True)

            # Define the output file path and save the DataFrame as a TSV file.
            output_file = 'conventional_combinations.tsv'
            primer_combinations_df.to_csv(output_file, sep='\t', index=False)
            print(f"Primer combinations saved to {output_file}")
        else:
            print("No valid primer combinations found.")

    except Exception as e:
        # If an error occurs during processing, print the error message.
        print(f"An error occurred: {e}")
        sys.exit(1)  # Exit the script with an error status.

if __name__ == "__main__":
    # Check for the correct number of command-line arguments when the script is executed.
    if len(sys.argv) != 2:
        # If the incorrect number of arguments is provided, display the usage message.
        print("Usage: python combine.py mapping_positions.tsv")
        sys.exit(1)  # Exit the script with an error status due to incorrect usage.
    else:
        # If the correct number of arguments is provided, proceed with processing.
        input_file = sys.argv[1]  # Assign the provided input file path to a variable.
        generate_primer_combinations(input_file)  # Call the function to generate primer combinations.
