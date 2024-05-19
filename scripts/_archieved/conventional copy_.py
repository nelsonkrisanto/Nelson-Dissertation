import pandas as pd
import sys
from tqdm import tqdm
import logging

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def generate_primer_combinations(input_file, regions_of_interest):
    try:
        logging.info(f"Loading mapping data from {input_file}")
        df_mapping = pd.read_csv(input_file, delimiter='\t')
        df_mapping['Primer'] = df_mapping['Primer'].str.strip().str.upper()
        logging.debug(f"Mapping DataFrame Head:\n{df_mapping.head()}")

        logging.info("Loading primer metadata from primer_metadata.tsv")
        df_metadata = pd.read_csv('primer_metadata.tsv', delimiter='\t')
        df_metadata['Primer'] = df_metadata['Primer'].str.strip().str.upper()
        logging.debug(f"Metadata DataFrame Head:\n{df_metadata.head()}")

        logging.info("Merging dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        logging.info(f"Merged dataframes successfully. Total primers: {len(merged_df)}")
        logging.debug(f"Merged DataFrame Head:\n{merged_df.head()}")

        fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')]
        rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')]

        primer_combinations = []

        logging.info("Searching for valid primer combinations")
        for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc="Processing Primers"):
            if pd.isna(fwd['Genotype']):
                logging.debug(f"Skipping forward primer {fwd['Primer']} due to missing genotype")
                continue

            matching_revs = rev_primers[((rev_primers['Reference'] == fwd['Reference']) & 
                                         ((rev_primers['Genotype'] == fwd['Genotype']) | 
                                          (rev_primers['Genotype'] == 'ALL'))) &
                                        (rev_primers['Start'] > fwd['Start'] + 300) & 
                                        (rev_primers['Start'] < fwd['Start'] + 1000)]
            
            logging.debug(f"Forward Primer: {fwd['Primer']} | Genotype: {fwd['Genotype']} | Matching Reverse Primers: {matching_revs['Primer'].tolist()}")

            if matching_revs.empty:
                logging.debug(f"No reverse primers found for forward primer: {fwd['Primer']} with genotype: {fwd['Genotype']} and reference: {fwd['Reference']}")
            
            for _, rev in matching_revs.iterrows():
                combination_name = f"{fwd['Primer']}_{rev['Primer']}"
                amplicon_length = rev['Start'] - fwd['Start']
                tm_max_diff = round(abs(fwd['Tm_max'] - rev['Tm_max']))
                tm_min_diff = round(abs(fwd['Tm_min'] - rev['Tm_min']))

                logging.debug(f"Evaluating combination: {combination_name} | Amplicon Length: {amplicon_length}, Tm Max Diff: {tm_max_diff}, Tm Min Diff: {tm_min_diff}")

                if tm_max_diff <= 20 and tm_min_diff <= 20:
                    for region, (start, end) in regions_of_interest.items():
                        if fwd['Start'] >= start and rev['End'] <= end:
                            logging.debug(f"Valid combination found: {combination_name} in region {region}")
                            primer_combinations.append({
                                'Forward_Primer': fwd['Primer'],
                                'Reverse_Primer': rev['Primer'],
                                'Combination_Name': combination_name,
                                'Reference': fwd['Reference'],
                                'Genotype': fwd['Genotype'],
                                'Region': region,
                                'Amplicon_Length': amplicon_length
                            })
                        else:
                            logging.debug(f"Region mismatch for {combination_name}: Fwd Start {fwd['Start']} | Rev End {rev['End']} | Region {region} ({start}-{end})")

        if primer_combinations:
            primer_combinations_df = pd.DataFrame(primer_combinations)
            primer_combinations_df.drop_duplicates(subset='Combination_Name', inplace=True)
            output_file = 'conventional_combinations.tsv'
            primer_combinations_df.to_csv(output_file, sep='\t', index=False)
            logging.info(f"Primer combinations saved to {output_file}")
        else:
            logging.info("No valid primer combinations found.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python combine.py mapping_positions.tsv")
        sys.exit(1)
    else:
        setup_logging()
        input_file = sys.argv[1]
        regions_of_interest = {
            'NS1': (1000, 2000),
            'NS3': (3000, 4000),
            'NS5': (5000, 7000)
        }
        generate_primer_combinations(input_file, regions_of_interest)
