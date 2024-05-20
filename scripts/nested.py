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

            if (fwd['Tm_min'] >= rev['Tm_min']) and (fwd['Tm_max'] <= rev['Tm_max']) and (tm_max_diff <= 5) and (tm_min_diff <= 5) and (gc_content_diff <= 10) and (fwd_homopolymer <= 4) and (rev_homopolymer <= 4):
                for region, (start, end) in regions_of_interest.items():
                    if fwd['Start'] >= start and rev['End'] <= end:
                        combinations.append({
                            'Forward_Primer': fwd['Primer'],
                            'Reverse_Primer': rev['Primer'],
                            'Combination_Name': combination_name,
                            'Reference': fwd['Reference'],
                            'Genotype': fwd['Genotype'],
                            'Amplicon_Length': rev['Amplicon_Length'],
                            'Round_Type': round_type,
                            'Region': region,
                            'First_Round_Start': fwd['Start'],
                            'First_Round_End': rev['End']
                        })
                        break
    return pd.DataFrame(combinations)

def match_inner_to_outer(outer_combos, inner_combos):
    matched_combinations = []
    for _, outer in outer_combos.iterrows():
        for _, inner in inner_combos.iterrows():
            if (inner['Reference'] == outer['Reference'] and 
                inner['Region'] == outer['Region'] and 
                inner['Genotype'] == outer['Genotype'] and 
                inner['First_Round_Start'] >= outer['First_Round_Start'] and 
                inner['First_Round_End'] <= outer['First_Round_End']):
                matched_combinations.append({
                    'Outer_Combination_Name': outer['Combination_Name'],
                    'Outer_Forward_Primer': outer['Forward_Primer'],
                    'Outer_Reverse_Primer': outer['Reverse_Primer'],
                    'Outer_Genotype': outer['Genotype'],
                    'Outer_Amplicon_Length': outer['Amplicon_Length'],
                    'Inner_Combination_Name': inner['Combination_Name'],
                    'Inner_Forward_Primer': inner['Forward_Primer'],
                    'Inner_Reverse_Primer': inner['Reverse_Primer'],
                    'Inner_Genotype': inner['Genotype'],
                    'Inner_Amplicon_Length': inner['Amplicon_Length'],
                    'Reference': outer['Reference'],
                    'Region': outer['Region']
                })
    return pd.DataFrame(matched_combinations)

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

        logging.info("Finding valid primer combinations for the second round")
        second_round_combinations = pd.DataFrame()
        for combo in first_round_combinations.itertuples(index=False):
            specific_fwd = fwd_primers[fwd_primers['Primer'] == combo.Forward_Primer].copy()
            specific_rev = rev_primers[rev_primers['Primer'] == combo.Reverse_Primer].copy()
            specific_fwd['First_Round_End'] = specific_rev['End'].values[0]
            temp_combos = find_combinations(specific_fwd, specific_rev, 0, inner_max_len, 0, regions_of_interest, 'second')
            
            # Ensure inner primers fall within the outer amplicon
            for _, inner_combo in temp_combos.iterrows():
                if inner_combo['First_Round_Start'] <= inner_combo['First_Round_End'] and inner_combo['First_Round_End'] >= inner_combo['First_Round_Start']:
                    second_round_combinations = pd.concat([second_round_combinations, pd.DataFrame([inner_combo])], ignore_index=True)

        all_combinations = pd.concat([first_round_combinations, second_round_combinations], ignore_index=True)
        if not all_combinations.empty:
            all_combinations.drop_duplicates(subset=['Combination_Name', 'Amplicon_Length', 'Reference'], inplace=True)
            
            all_combinations_file = 'all_nested_primer_combinations.tsv'
            all_combinations.to_csv(all_combinations_file, sep='\t', index=False)
            logging.info(f"All nested primer combinations saved to {all_combinations_file}")
            
            top_200_combinations_file = 'top_200_nested_primer_combinations.tsv'
            top_200_combinations_df = all_combinations.head(200)
            top_200_combinations_df.to_csv(top_200_combinations_file, sep='\t', index=False)
            logging.info(f"Top 200 nested primer combinations saved to {top_200_combinations_file}")

            # Match inner combinations to outer combinations
            logging.info("Matching inner combinations to outer combinations")
            matched_combinations = match_inner_to_outer(first_round_combinations, second_round_combinations)
            matched_combinations_file = 'matched_nested_primer_combinations.tsv'
            matched_combinations.to_csv(matched_combinations_file, sep='\t', index=False)
            logging.info(f"Matched nested primer combinations saved to {matched_combinations_file}")

        else:
            logging.info("No valid nested primer combinations found.")
    except Exception as e:
                logging.error(f"An error occurred: {e}")
                sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py mapping_positions.tsv primer_metadata.tsv outer_min_len inner_max_len")
        sys.exit(1)
    else:
        setup_logging()
        mapping_file, metadata_file, outer_min_len, inner_max_len = sys.argv[1:5]
        generate_nested_primer_combinations(mapping_file, metadata_file, int(outer_min_len), int(inner_max_len))
