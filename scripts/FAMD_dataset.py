import pandas as pd
import sys
from tqdm import tqdm
import logging
from itertools import groupby

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_gc_content(sequence):
    return round((sequence.count('G') + sequence.count('C')) / len(sequence) * 100, 2)

def check_primer_conditions(primer):
    length = len(primer['Sequence'])
    gc_content = calculate_gc_content(primer['Sequence'])
    homopolymer_run = max([len(list(g)) for k, g in groupby(primer['Sequence'])])
    if length < 18 or length > 25:
        return False
    if gc_content < 40 or gc_content > 60:
        return False
    if homopolymer_run > 4:
        return False
    return True

# Function to extract individual primer details
def extract_primer_details(mapping_file, metadata_file, combinations_file, output_file, regions_of_interest):
    try:
        # Load the mapping data from the provided TSV file
        logging.info(f"Loading mapping data from {mapping_file}")
        df_mapping = pd.read_csv(mapping_file, delimiter='\t')
        df_mapping['Primer'] = df_mapping['Primer'].str.strip().str.upper()
        logging.debug(f"Mapping DataFrame Head:\n{df_mapping.head()}")

        # Load the primer metadata from the provided TSV file
        logging.info(f"Loading primer metadata from {metadata_file}")
        df_metadata = pd.read_csv(metadata_file, delimiter='\t')
        df_metadata['Primer'] = df_metadata['Primer'].str.strip().str.upper()
        df_metadata['GC_Content'] = df_metadata['Sequence'].apply(calculate_gc_content)
        logging.debug(f"Metadata DataFrame Head:\n{df_metadata.head()}")

        # Load the primer combinations from the provided TSV file
        logging.info(f"Loading primer combinations from {combinations_file}")
        df_combinations = pd.read_csv(combinations_file, delimiter='\t', header=None, names=['Sequence_F', 'Sequence_R', 'Primer'])
        df_combinations['Primer'] = df_combinations['Primer'].str.strip().str.upper()
        logging.debug(f"Combinations DataFrame Head:\n{df_combinations.head()}")

        # Merge the mapping and metadata dataframes based on the 'Primer' column
        logging.info("Merging mapping and metadata dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        logging.info(f"Merged dataframes successfully. Total primers: {len(merged_df)}")
        logging.debug(f"Merged DataFrame Head:\n{merged_df.head()}")

        primer_details = []

        # Start extracting primer details
        logging.info("Extracting primer details")
        for _, primer in tqdm(merged_df.iterrows(), total=merged_df.shape[0], desc="Processing Primers"):
            if pd.isna(primer['Genotype_x']) or not check_primer_conditions(primer):
                logging.debug(f"Skipping primer {primer['Primer']} due to missing genotype or failing conditions")
                continue

            orientation = 'Forward' if primer['Primer'].endswith('_F') else 'Reverse'

            for region, (start, end) in regions_of_interest.items():
                if primer['Start'] >= start and primer['End'] <= end:
                    logging.debug(f"Valid primer found: {primer['Primer']} in region {region}")
                    primer_details.append({
                        'Primer_Name': primer['Primer'],
                        'Sequence': primer['Sequence'],
                        'Genotype': primer['Genotype_x'],
                        'Orientation': orientation,
                        'Start': primer['Start'],
                        'End': primer['End'],
                        'Amplicon_Length': primer['End'] - primer['Start'],
                        'Region': region,
                        'GC_Content': primer['GC_Content'],
                        'Tm_min': primer['Tm_min'],
                        'Tm_max': primer['Tm_max']
                    })

        # Save the extracted primer details if any are found
        if primer_details:
            primer_details_df = pd.DataFrame(primer_details)
            primer_details_df.drop_duplicates(subset=['Primer_Name'], inplace=True)
            
            # Rearrange columns as requested
            primer_details_df = primer_details_df[[
                'Primer_Name', 'Sequence', 'Genotype', 'Orientation', 
                'Start', 'End', 'Amplicon_Length', 'Region', 
                'GC_Content', 'Tm_min', 'Tm_max'
            ]]

            # Save the primer details to a file
            primer_details_df.to_csv(output_file, sep='\t', index=False)
            
            logging.info(f"Individual primer details saved to {output_file}")
        else:
            logging.info("No valid primers found.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)

# Main entry point of the script
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python analyze_primers.py mapping_file metadata_file combinations_file output_file")
        sys.exit(1)
    else:
        setup_logging()
        mapping_file = sys.argv[1]
        metadata_file = sys.argv[2]
        combinations_file = sys.argv[3]
        output_file = sys.argv[4]
        regions_of_interest = {
            'NS1': (1000, 2000),
            'NS3': (3000, 4000),
            'NS5': (5000, 7000)
        }
        extract_primer_details(mapping_file, metadata_file, combinations_file, output_file, regions_of_interest)
