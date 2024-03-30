import pandas as pd
from tqdm import tqdm
import sys

def prefilter_data(merged_df):
    """
    Separates the merged DataFrame into two DataFrames containing forward and reverse primers
    """
    # Filters out forward primers and creates a copy to avoid SettingWithCopyWarning.
    fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')].copy()
    # Filters out reverse primers and creates a copy.
    rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')].copy()
    return fwd_primers, rev_primers

def find_semi_nested_combinations(fwd_primers, rev_primers, is_forward_reused, min_length, max_length):
    """
    Identifies valid primer combinations for semi-nested PCR, distinguishing between scenarios
    where a forward primer is reused (and a new reverse primer is used) and vice versa.
    """
    # Determines which set of primers are reused and which are new based on whether a forward primer is being reused.
    reused_primers, new_primers = (fwd_primers, rev_primers) if is_forward_reused else (rev_primers, fwd_primers)

    combinations = []  # Initializes a list to store valid primer combinations.
    # Iterates over reused primers with progress tracking.
    for _, reused_primer in tqdm(reused_primers.iterrows(), total=reused_primers.shape[0], desc="Processing Reused Primers"):
        # Filters new primers by matching references with the reused primer.
        matching_new_primers = new_primers[new_primers['Reference'] == reused_primer['Reference']]
        
        # Iterates over matching new primers to check for valid combinations.
        for _, new_primer in matching_new_primers.iterrows():
            # Calculates the amplicon length for potential semi-nested PCR combinations.
            amplicon_length = abs(new_primer['Start'] - reused_primer['Start'])
            # Checks if the amplicon length falls within the specified range.
            if min_length <= amplicon_length <= max_length:
                # If valid, adds the combination to the list with relevant details.
                combinations.append({
                    'Reused_Primer': reused_primer['Primer'],
                    'New_Primer': new_primer['Primer'],
                    'Combination_Name': f"{reused_primer['Primer']}_{new_primer['Primer']}",
                    'Reference': reused_primer['Reference'],
                    'Amplicon_Length': amplicon_length
                })
    
    return combinations  # Returns the list of valid primer combinations.

def generate_semi_nested_primer_combinations(mapping_positions_path, primer_metadata_path):
    """
    Main function that generate the semi-nested PCR primer combinations.
    It processes both forward and reverse primers to find suitable pairs for semi-nested PCR.
    """
    try:
        # Loads mapping positions and primer metadata into DataFrames.
        df_mapping = pd.read_csv(mapping_positions_path, sep='\t')
        df_metadata = pd.read_csv(primer_metadata_path, sep='\t')

        # Merges the mapping and metadata DataFrames on the 'Primer' column.
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        # Pre-filters the merged DataFrame to separate forward and reverse primers.
        fwd_primers, rev_primers = prefilter_data(merged_df)

        # Finds valid semi-nested combinations where forward primers are reused.
        forward_combinations = find_semi_nested_combinations(fwd_primers, rev_primers, True, 100, 500)
        print(f"Found {len(forward_combinations)} valid semi-nested combinations with forward primers reused.")

        # Finds valid semi-nested combinations where reverse primers are reused.
        reverse_combinations = find_semi_nested_combinations(fwd_primers, rev_primers, False, 100, 500)
        print(f"Found {len(reverse_combinations)} valid semi-nested combinations with reverse primers reused.")

        # Processes the valid combinations for forward reused primers, if any.
        if forward_combinations:
            forward_df = pd.DataFrame(forward_combinations)
            forward_df.drop_duplicates(subset=['Combination_Name', 'Amplicon_Length', 'Reference'], inplace=True)
            forward_file = 'semi_nested_combinations_forward.tsv'
            forward_df.to_csv(forward_file, sep='\t', index=False)
            print(f"Semi-nested primer combinations with forward primers reused saved to {forward_file}")

        # Processes the valid combinations for reverse reused primers, if any.
        if reverse_combinations:
            reverse_df = pd.DataFrame(reverse_combinations)
            reverse_df.drop_duplicates(subset=['Combination_Name', 'Amplicon_Length', 'Reference'], inplace=True)
            reverse_file = 'semi_nested_combinations_reverse.tsv'
            reverse_df.to_csv(reverse_file, sep='\t', index=False)
            print(f"Semi-nested primer combinations with reverse primers reused saved to {reverse_file}")

        # If no valid combinations are found for either case, prints a message.
        if not forward_combinations and not reverse_combinations:
            print("No valid semi-nested primer combinations found.")

    # Handles various exceptions that may arise during execution.
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
    # Checks for the correct number of command-line arguments.
    if len(sys.argv) != 3:
        print("Usage: python script.py mapping_positions.tsv primer_metadata.tsv")
        sys.exit(1)
    else:
        # Assigns file paths based on command-line arguments.
        mapping_positions_path = sys.argv[1]
        primer_metadata_path = sys.argv[2]
        # Calls the main function to generate semi-nested PCR primer combinations.
        generate_semi_nested_primer_combinations(mapping_positions_path, primer_metadata_path)
