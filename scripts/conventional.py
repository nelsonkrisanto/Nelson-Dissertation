import pandas as pd
import sys
from tqdm import tqdm
import logging

# Function to set up logging configuration
def setup_logging():
    # Configure logging to display the time, log level, and message
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Function to generate primer combinations
def generate_primer_combinations(input_file, regions_of_interest):
    try:
        # Load the mapping data from the provided TSV file
        logging.info(f"Loading mapping data from {input_file}")
        df_mapping = pd.read_csv(input_file, delimiter='\t')
        # Clean and standardize the 'Primer' column by stripping whitespace and converting to uppercase
        df_mapping['Primer'] = df_mapping['Primer'].str.strip().str.upper()
        logging.debug(f"Mapping DataFrame Head:\n{df_mapping.head()}")

        # Load the primer metadata from 'primer_metadata.tsv'
        logging.info("Loading primer metadata from primer_metadata.tsv")
        df_metadata = pd.read_csv('primer_metadata.tsv', delimiter='\t')
        # Clean and standardize the 'Primer' column by stripping whitespace and converting to uppercase
        df_metadata['Primer'] = df_metadata['Primer'].str.strip().str.upper()
        logging.debug(f"Metadata DataFrame Head:\n{df_metadata.head()}")

        # Merge the mapping and metadata dataframes based on the 'Primer' column
        logging.info("Merging dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        logging.info(f"Merged dataframes successfully. Total primers: {len(merged_df)}")
        logging.debug(f"Merged DataFrame Head:\n{merged_df.head()}")

        # Separate the merged dataframe into forward and reverse primers
        fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')]
        rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')]

        primer_combinations = []

        # Start searching for valid primer combinations
        logging.info("Searching for valid primer combinations")
        for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc="Processing Primers"):
            # Skip forward primers without a genotype
            if pd.isna(fwd['Genotype_x']):
                logging.debug(f"Skipping forward primer {fwd['Primer']} due to missing genotype")
                continue

            # Find matching reverse primers based on specific criteria
            matching_revs = rev_primers[((rev_primers['Reference'] == fwd['Reference']) & 
                                         ((rev_primers['Genotype_x'] == fwd['Genotype_x']) | 
                                          (rev_primers['Genotype_x'] == 'ALL'))) &
                                        (rev_primers['Start'] > fwd['Start'] + 300) & 
                                        (rev_primers['Start'] < fwd['Start'] + 1050)]  # Increased distance to account for the 50bp difference
            
            logging.debug(f"Forward Primer: {fwd['Primer']} | Genotype: {fwd['Genotype_x']} | Matching Reverse Primers: {matching_revs['Primer'].tolist()}")

            # Log if no matching reverse primers are found
            if matching_revs.empty:
                logging.debug(f"No reverse primers found for forward primer: {fwd['Primer']} with genotype: {fwd['Genotype_x']} and reference: {fwd['Reference']}")
            
            # Evaluate each matching reverse primer
            for _, rev in matching_revs.iterrows():
                combination_name = f"{fwd['Primer']}_{rev['Primer']}"
                amplicon_length = rev['Start'] - fwd['Start']
                tm_max_diff = abs(fwd['Tm_max'] - rev['Tm_max'])
                tm_min_diff = abs(fwd['Tm_min'] - rev['Tm_min'])

                logging.debug(f"Evaluating combination: {combination_name} | Amplicon Length: {amplicon_length}, Tm Max Diff: {tm_max_diff}, Tm Min Diff: {tm_min_diff}")

                # Check if the combination meets the Tm and region criteria
                if (fwd['Tm_min'] >= rev['Tm_min']) and (fwd['Tm_max'] <= rev['Tm_max']) and (tm_max_diff <= 5) and (tm_min_diff <= 5):
                    for region, (start, end) in regions_of_interest.items():
                        if fwd['Start'] >= start and rev['End'] <= end:
                            logging.debug(f"Valid combination found: {combination_name} in region {region}")
                            primer_combinations.append({
                                'Forward_Primer': fwd['Primer'],
                                'Reverse_Primer': rev['Primer'],
                                'Combination_Name': combination_name,
                                'Reference': fwd['Reference'],
                                'Genotype': fwd['Genotype_x'],
                                'Region': region,
                                'Amplicon_Length': amplicon_length
                            })
                        else:
                            logging.debug(f"Region mismatch for {combination_name}: Fwd Start {fwd['Start']} | Rev End {rev['End']} | Region {region} ({start}-{end})")

        # Save the valid primer combinations if any are found
        if primer_combinations:
            primer_combinations_df = pd.DataFrame(primer_combinations)
            primer_combinations_df.drop_duplicates(subset=['Combination_Name', 'Amplicon_Length', 'Reference'], inplace=True)
            
            # Save all primer combinations to a file
            all_combinations_file = 'all_conv_primer_combinations.tsv'
            primer_combinations_df.to_csv(all_combinations_file, sep='\t', index=False)
            logging.info(f"All primer combinations saved to {all_combinations_file}")
            
            # Save the top 200 primer combinations to a separate file
            top_200_combinations_file = 'top_200_conv_primer_combinations.tsv'
            top_200_combinations_df = primer_combinations_df.head(200)
            top_200_combinations_df.to_csv(top_200_combinations_file, sep='\t', index=False)
            logging.info(f"Top 200 primer combinations saved to {top_200_combinations_file}")
        else:
            logging.info("No valid primer combinations found.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)

# Main entry point of the script
if __name__ == "__main__":
    # Check if the script is called with the correct number of arguments
    if len(sys.argv) != 2:
        print("Usage: python combine.py mapping_positions.tsv")
        sys.exit(1)
    else:
        # Set up logging
        setup_logging()
        # Get the input file path from the command line arguments
        input_file = sys.argv[1]
        # Define the regions of interest
        regions_of_interest = {
            'NS1': (1000, 2000),
            'NS3': (3000, 4000),
            'NS5': (5000, 7000)
        }
        # Generate primer combinations
        generate_primer_combinations(input_file, regions_of_interest)
