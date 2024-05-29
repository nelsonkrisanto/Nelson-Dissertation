import pandas as pd
from tqdm import tqdm
import sys
import logging

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def prefilter_data(merged_df):
    fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')].copy()
    rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')].copy()
    return fwd_primers, rev_primers

def check_primer_conditions(primer):
    length = len(primer['Sequence'])
    gc_content = (primer['Sequence'].count('G') + primer['Sequence'].count('C')) / length * 100
    if length < 18 or length > 25:
        return False
    if gc_content < 40 or gc_content > 60:
        return False
    if 'AAAA' in primer['Sequence'] or 'TTTT' in primer['Sequence'] or 'GGGG' in primer['Sequence'] or 'CCCC' in primer['Sequence']:
        return False
    return True

def find_combinations(fwd_primers, rev_primers, min_length, max_length, min_primer_distance, regions_of_interest, round_type='first'):
    combinations = []
    for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc=f"Processing {round_type.capitalize()} Round Forward Primers"):
        if pd.isna(fwd['Genotype']) or not check_primer_conditions(fwd):
            logging.debug(f"Skipping forward primer {fwd['Primer']} due to missing genotype or failing conditions")
            continue

        matching_revs = rev_primers[(rev_primers['Reference'] == fwd['Reference']) & 
                                    ((rev_primers['Genotype'] == fwd['Genotype']) | 
                                     (rev_primers['Genotype'] == 'ALL')) &
                                    (rev_primers['Start'] > fwd['Start'] + min_primer_distance)]
        matching_revs['Amplicon_Length'] = matching_revs['Start'] - fwd['Start']
        
        valid_revs = matching_revs[(matching_revs['Amplicon_Length'] >= min_length) & 
                                   (matching_revs['Amplicon_Length'] <= max_length)]
        
        for _, rev in valid_revs.iterrows():
            if not check_primer_conditions(rev):
                continue

            combination_name = f"{fwd['Primer']}_{rev['Primer']}"
            tm_max_diff = abs(fwd['Tm_max'] - rev['Tm_max'])
            tm_min_diff = abs(fwd['Tm_min'] - rev['Tm_min'])
            
            if (fwd['Tm_min'] >= rev['Tm_min']) and (fwd['Tm_max'] <= rev['Tm_max']) and (tm_max_diff <= 3) and (tm_min_diff <= 3):
                for region, (start, end) in regions_of_interest.items():
                    if fwd['Start'] >= start and rev['End'] <= end:
                        combinations.append({
                            'Forward_Primer': fwd['Sequence'],
                            'Reverse_Primer': rev['Sequence'],
                            'Combination_Name': combination_name.replace(' ', '_'),
                            'Round_Type': round_type,
                            'First_Round_Start': fwd['Start'],
                            'First_Round_End': rev['End'],
                            'Reference': fwd['Reference']
                        })
                        break
    return pd.DataFrame(combinations)

def generate_semi_nested_primer_combinations(mapping_file, metadata_file, outer_min_len, inner_max_len, regions_of_interest):
    try:
        logging.info(f"Loading mapping data from {mapping_file}")
        df_mapping = pd.read_csv(mapping_file, delimiter='\t')
        df_mapping['Primer'] = df_mapping['Primer'].str.strip().str.upper()
        logging.debug(f"Mapping DataFrame Head:\n{df_mapping.head()}")

        logging.info(f"Loading primer metadata from {metadata_file}")
        df_metadata = pd.read_csv(metadata_file, delimiter='\t')
        df_metadata['Primer'] = df_metadata['Primer'].str.strip().str.upper()
        df_metadata['GC_Content'] = df_metadata['Sequence'].apply(lambda seq: (seq.count('G') + seq.count('C')) / len(seq) * 100)
        logging.debug(f"Metadata DataFrame Head:\n{df_metadata.head()}")

        logging.info("Merging dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        logging.info(f"Merged dataframes successfully. Total primers: {len(merged_df)}")
        logging.debug(f"Merged DataFrame Head:\n{merged_df.head()}")

        # Rename columns to avoid issues with suffixes
        merged_df.rename(columns={'Genotype_x': 'Genotype'}, inplace=True)

        fwd_primers, rev_primers = prefilter_data(merged_df)

        min_primer_distance = 50  # Minimum distance between primer start positions for both rounds

        logging.info("Finding valid primer combinations for the first round")
        first_round_combinations = find_combinations(fwd_primers, rev_primers, outer_min_len, float('inf'), min_primer_distance, regions_of_interest)

        logging.info("Finding valid primer combinations for the second round")
        second_round_combinations = pd.DataFrame()
        for combo in first_round_combinations.itertuples(index=False):
            specific_fwd = fwd_primers[(fwd_primers['Reference'] == combo.Reference) & (fwd_primers['Start'] > combo.First_Round_Start) & (fwd_primers['Start'] < combo.First_Round_End)].copy()
            specific_rev = rev_primers[(rev_primers['Reference'] == combo.Reference) & (rev_primers['End'] < combo.First_Round_End) & (rev_primers['End'] > combo.First_Round_Start)].copy()
            if specific_fwd.empty or specific_rev.empty:
                logging.debug(f"No matching forward or reverse primers for second round: {combo.Forward_Primer}, {combo.Reverse_Primer}")
                continue
            
            # Find new forward or reverse primer for semi-nested PCR
            new_fwd_combinations = find_combinations(specific_fwd, rev_primers, 0, inner_max_len, 0, regions_of_interest, 'second')
            new_rev_combinations = find_combinations(fwd_primers, specific_rev, 0, inner_max_len, 0, regions_of_interest, 'second')
            
            second_round_combinations = pd.concat([second_round_combinations, new_fwd_combinations, new_rev_combinations], ignore_index=True)

        if not first_round_combinations.empty:
            first_round_combinations.drop_duplicates(subset=['Combination_Name'], inplace=True)
            outer_combinations_file = 'semi_nested_outer_primer_combinations.tsv'
            with open(outer_combinations_file, 'w') as f:
                for i, row in first_round_combinations.iterrows():
                    f.write(f"{row['Forward_Primer']}\t{row['Reverse_Primer']}\t{row['Combination_Name']}\n")
            logging.info(f"Formatted outer primer combinations saved to {outer_combinations_file}")

        if not second_round_combinations.empty:
            second_round_combinations.drop_duplicates(subset=['Combination_Name'], inplace=True)
            inner_combinations_file = 'semi_nested_inner_primer_combinations.tsv'
            with open(inner_combinations_file, 'w') as f:
                for i, row in second_round_combinations.iterrows():
                    f.write(f"{row['Forward_Primer']}\t{row['Reverse_Primer']}\t{row['Combination_Name']}\n")
            logging.info(f"Formatted inner primer combinations saved to {inner_combinations_file}")
        else:
            logging.info("No valid inner primer combinations found.")
        
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
        regions_of_interest = {
            'NS1': (1000, 2000),
            'NS3': (3000, 4000),
            'NS5': (5000, 7000)
        }
        generate_semi_nested_primer_combinations(mapping_file, metadata_file, int(outer_min_len), int(inner_max_len), regions_of_interest)
