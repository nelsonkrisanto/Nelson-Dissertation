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

# Function to generate primer combinations
def generate_primer_combinations(mapping_file, metadata_file, output_file, regions_of_interest):
    try:
        setup_logging()

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

        # Merge the mapping and metadata dataframes based on the 'Primer' column
        logging.info("Merging dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        logging.info(f"Merged dataframes successfully. Total primers: {len(merged_df)}")
        logging.debug(f"Merged DataFrame Head:\n{merged_df.head()}")

        # Rename columns to avoid issues with suffixes
        merged_df.rename(columns={'Genotype_x': 'Genotype', 'Reference_x': 'Reference'}, inplace=True)

        # Separate the merged dataframe into forward and reverse primers
        fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')]
        rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')]

        primer_combinations = []

        # Start searching for valid primer combinations
        logging.info("Searching for valid primer combinations")
        for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc="Processing Primers"):
            if pd.isna(fwd['Genotype']) or not check_primer_conditions(fwd):
                logging.debug(f"Skipping forward primer {fwd['Primer']} due to missing genotype or failing conditions")
                continue

            matching_revs = rev_primers[(rev_primers['Reference'] == fwd['Reference']) & 
                                        ((rev_primers['Genotype'] == fwd['Genotype']) | 
                                         (rev_primers['Genotype'] == 'ALL')) &
                                        (rev_primers['Start'] > fwd['Start'] + 300) & 
                                        (rev_primers['Start'] < fwd['Start'] + 1200)]
            matching_revs['Amplicon_Length'] = matching_revs['Start'] - fwd['Start']

            logging.debug(f"Forward Primer: {fwd['Primer']} | Genotype: {fwd['Genotype']} | Matching Reverse Primers: {matching_revs['Primer'].tolist()}")

            if matching_revs.empty:
                logging.debug(f"No reverse primers found for forward primer: {fwd['Primer']} with genotype: {fwd['Genotype']} and reference: {fwd['Reference']}")

            for _, rev in matching_revs.iterrows():
                if not check_primer_conditions(rev):
                    continue

                combination_name = f"{fwd['Primer']}_{rev['Primer']}"
                tm_max_diff = abs(fwd['Tm_max'] - rev['Tm_max'])
                tm_min_diff = abs(fwd['Tm_min'] - rev['Tm_min'])
                gc_content_diff = abs(fwd['GC_Content'] - rev['GC_Content'])

                logging.debug(f"Evaluating combination: {combination_name} | Amplicon Length: {rev['Amplicon_Length']}, Tm Max Diff: {tm_max_diff}, Tm Min Diff: {tm_min_diff}, GC Content Diff: {gc_content_diff}")

                if (fwd['Tm_min'] >= rev['Tm_min']) and (fwd['Tm_max'] <= rev['Tm_max']) and (tm_max_diff <= 3) and (tm_min_diff <= 3) and (gc_content_diff <= 10):
                    for region, (start, end) in regions_of_interest.items():
                        if fwd['Start'] >= start and rev['End'] <= end:
                            logging.debug(f"Valid combination found: {combination_name} in region {region}")
                            primer_combinations.append({
                                'Amplicon_Length': rev['Amplicon_Length'],
                                'Combination_Name': combination_name.replace(' ', '_'),
                                'End': rev['End'],
                                'Forward_Primer': fwd['Primer'],
                                'GC_Content_Forward': fwd['GC_Content'],
                                'GC_Content_Reverse': rev['GC_Content'],
                                'Genotype_Forward': fwd['Genotype'],
                                'Genotype_Reverse': rev['Genotype'],
                                'Reference': fwd['Reference'],
                                'Region': region,
                                'Reverse_Primer': rev['Primer'],
                                'Sequence_Forward': fwd['Sequence'],
                                'Sequence_Reverse': rev['Sequence'],
                                'Start': fwd['Start'],
                                'Tm_max_Forward': fwd['Tm_max'],
                                'Tm_max_Reverse': rev['Tm_max'],
                                'Tm_min_Forward': fwd['Tm_min'],
                                'Tm_min_Reverse': rev['Tm_min']
                            })

        # Save the valid primer combinations if any are found
        if primer_combinations:
            primer_combinations_df = pd.DataFrame(primer_combinations)
            primer_combinations_df.drop_duplicates(subset=['Combination_Name'], inplace=True)
            
            # Save the primer combinations to a file
            primer_combinations_df.to_csv(output_file, sep='\t', index=False)
            
            logging.info(f"Formatted conventional primer combinations saved to {output_file}")
        else:
            logging.info("No valid primer combinations found.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)

# Main entry point of the script
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python combine.py mapping_positions.tsv primer_metadata.tsv output_file.tsv")
        sys.exit(1)
    else:
        setup_logging()
        mapping_file = sys.argv[1]
        metadata_file = sys.argv[2]
        output_file = sys.argv[3]
        regions_of_interest = {
            'NS1': (1000, 2000),
            'NS3': (3000, 4000),
            'NS5': (5000, 7000)
        }
        generate_primer_combinations(mapping_file, metadata_file, output_file, regions_of_interest)
