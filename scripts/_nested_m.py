import pandas as pd
from tqdm import tqdm
import sys
import logging
from itertools import groupby

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def prefilter_data(merged_df):
    fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')].copy()
    rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')].copy()
    return fwd_primers, rev_primers

def calculate_gc_content(sequence):
    return round((sequence.count('G') + sequence.count('C')) / len(sequence) * 100, 2)

def find_combinations(fwd_primers, rev_primers, min_length, max_length, min_primer_distance, regions_of_interest, round_type='first'):
    combinations = []
    for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc=f"Processing {round_type.capitalize()} Round Forward Primers"):
        if pd.isna(fwd['Genotype']):
            logging.debug(f"Skipping forward primer {fwd['Primer']} due to missing genotype")
            continue

        matching_revs = rev_primers[(rev_primers['Reference'] == fwd['Reference']) & 
                                    ((rev_primers['Genotype'] == fwd['Genotype']) | 
                                     (rev_primers['Genotype'] == 'ALL')) &
                                    (rev_primers['Start'] > fwd['Start'] + min_primer_distance)]
        matching_revs['Amplicon_Length'] = matching_revs['Start'] - fwd['Start']
        
        valid_revs = matching_revs[(matching_revs['Amplicon_Length'] >= min_length) & 
                                   (matching_revs['Amplicon_Length'] <= max_length)]
        
        for _, rev in valid_revs.iterrows():
            combination_name = f"{fwd['Primer']}_{rev['Primer']}"
            tm_max_diff = abs(fwd['Tm_max'] - rev['Tm_max'])
            tm_min_diff = abs(fwd['Tm_min'] - rev['Tm_min'])
            gc_content_diff = abs(fwd['GC_Content'] - rev['GC_Content'])
            fwd_homopolymer = max([len(list(g)) for k, g in groupby(fwd['Sequence'])])
            rev_homopolymer = max([len(list(g)) for k, g in groupby(rev['Sequence'])])

            if (fwd['Tm_min'] >= rev['Tm_min']) and (fwd['Tm_max'] <= rev['Tm_max']) and (tm_max_diff <= 3) and (tm_min_diff <= 3) and (gc_content_diff <= 10) and (fwd_homopolymer <= 4) and (rev_homopolymer <= 4):
                for region, (start, end) in regions_of_interest.items():
                    if fwd['Start'] >= start and rev['End'] <= end:
                        combinations.append({
                            'Forward_Primer': fwd['Sequence'],
                            'Reverse_Primer': rev['Sequence'],
                            'Combination_Name': combination_name.replace(' ', '_'),
                        })
                        break
    return pd.DataFrame(combinations)

def generate_nested_primer_combinations(mapping_file, metadata_file, outer_min_len, inner_max_len):
    try:
        logging.info(f"Loading mapping data from {mapping_file}")
        df_mapping = pd.read_csv(mapping_file, delimiter='\t')
        df_mapping['Primer'] = df_mapping['Primer'].str.strip().str.upper()
        logging.debug(f"Mapping DataFrame Head:\n{df_mapping.head()}")

        logging.info(f"Loading primer metadata from {metadata_file}")
        df_metadata = pd.read_csv(metadata_file, delimiter='\t')
        df_metadata['Primer'] = df_metadata['Primer'].str.strip().str.upper()
        df_metadata['GC_Content'] = df_metadata['Sequence'].apply(calculate_gc_content)
        logging.debug(f"Metadata DataFrame Head:\n{df_metadata.head()}")

        logging.info("Merging dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        logging.info(f"Merged dataframes successfully. Total primers: {len(merged_df)}")
        logging.debug(f"Merged DataFrame Head:\n{merged_df.head()}")

        # Rename columns to avoid issues with suffixes
        merged_df.rename(columns={'Genotype_x': 'Genotype'}, inplace=True)

        fwd_primers, rev_primers = prefilter_data(merged_df)

        min_primer_distance = 50  # Minimum distance between primer start positions for both rounds

        regions_of_interest = {
            'NS1': (1000, 2000),
            'NS3': (3000, 4000),
            'NS5': (5000, 7000)
        }

        logging.info("Finding valid primer combinations for the first round")
        first_round_combinations = find_combinations(fwd_primers, rev_primers, outer_min_len, float('inf'), min_primer_distance, regions_of_interest)

        if first_round_combinations.empty:
            logging.info("No valid first-round primer combinations found.")
            return

        logging.info("Finding valid primer combinations for the second round")
        second_round_combinations = pd.DataFrame()
        for combo in first_round_combinations.itertuples(index=False):
            specific_fwd = fwd_primers[fwd_primers['Primer'] == combo.Forward_Primer].copy()
            specific_rev = rev_primers[rev_primers['Primer'] == combo.Reverse_Primer].copy()
            if specific_rev.empty:
                continue
            specific_fwd['First_Round_End'] = specific_rev['End'].values[0]
            temp_combos = find_combinations(specific_fwd, specific_rev, 0, inner_max_len, 0, regions_of_interest, 'second')
            
            # Ensure inner primers fall within the outer amplicon
            for _, inner_combo in temp_combos.iterrows():
                if inner_combo['First_Round_Start'] <= inner_combo['First_Round_End'] and inner_combo['First_Round_End'] >= inner_combo['First_Round_Start']:
                    second_round_combinations = pd.concat([second_round_combinations, pd.DataFrame([inner_combo])], ignore_index=True)

        all_combinations = pd.concat([first_round_combinations, second_round_combinations], ignore_index=True)
        if not all_combinations.empty:
            all_combinations.drop_duplicates(subset=['Combination_Name'], inplace=True)
            
            # Save primer combinations to a file in the specified format
            output_file = 'nested_primer_combinations_m.tsv'
            with open(output_file, 'w') as f:
                for i, row in all_combinations.iterrows():
                    f.write(f"{row['Forward_Primer']}\t{row['Reverse_Primer']}\t{row['Combination_Name']}\n")
            logging.info(f"Formatted nested primer combinations saved to {output_file}")

        else:
            logging.info("No valid primer combinations found.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)

# Main entry point of the script
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py mapping_positions.tsv primer_metadata.tsv outer_min_len inner_max_len")
        sys.exit(1)
    else:
        setup_logging()
        mapping_file, metadata_file, outer_min_len, inner_max_len = sys.argv[1:5]
        generate_nested_primer_combinations(mapping_file, metadata_file, int(outer_min_len), int(inner_max_len))
