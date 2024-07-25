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

def generate_semi_nested_primer_combinations(input_file, regions_of_interest):
    try:
        logging.info(f"Loading mapping data from {input_file}")
        df_mapping = pd.read_csv(input_file, delimiter='\t')
        df_mapping['Primer'] = df_mapping['Primer'].str.strip().str.upper()
        logging.debug(f"Mapping DataFrame Head:\n{df_mapping.head()}")

        logging.info("Loading primer metadata from primer_metadata.tsv")
        df_metadata = pd.read_csv('primer_metadata.tsv', delimiter='\t')
        df_metadata['Primer'] = df_metadata['Primer'].str.strip().str.upper()
        df_metadata['GC_Content'] = df_metadata['Sequence'].apply(calculate_gc_content)
        logging.debug(f"Metadata DataFrame Head:\n{df_metadata.head()}")

        logging.info("Merging dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        logging.info(f"Merged dataframes successfully. Total primers: {len(merged_df)}")
        logging.debug(f"Merged DataFrame Head:\n{merged_df.head()}")

        merged_df.rename(columns={'Genotype_x': 'Genotype'}, inplace=True)

        fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')]
        rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')]

        outer_primer_combinations = []
        semi_nested_primer_combinations = []

        logging.info("Searching for valid outer primer combinations")
        for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc="Processing Outer Primers"):
            if pd.isna(fwd['Genotype']) or not check_primer_conditions(fwd):
                logging.debug(f"Skipping forward primer {fwd['Primer']} due to missing genotype or failing conditions")
                continue

            matching_revs = rev_primers[(rev_primers['Reference'] == fwd['Reference']) & 
                                        ((rev_primers['Genotype'] == fwd['Genotype']) | 
                                         (rev_primers['Genotype'] == 'ALL')) &
                                        (rev_primers['Start'] > fwd['Start'] + 500) & 
                                        (rev_primers['Start'] < fwd['Start'] + 2000)]
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

                if (tm_max_diff <= 3) and (tm_min_diff <= 3) and (gc_content_diff <= 10):
                    for region, (start, end) in regions_of_interest.items():
                        if fwd['Start'] >= start and rev['End'] <= end:
                            logging.debug(f"Valid combination found: {combination_name} in region {region}")
                            outer_primer_combinations.append({
                                'Forward_Primer': fwd['Primer'],
                                'Forward_Sequence': fwd['Sequence'],
                                'Reverse_Primer': rev['Primer'],
                                'Reverse_Sequence': rev['Sequence'],
                                'Combination_Name': combination_name.replace(' ', '_'),
                                'Forward_Start': fwd['Start'],
                                'Reverse_End': rev['Start'],
                                'Tm_max': fwd['Tm_max'],
                                'Tm_min': fwd['Tm_min'],
                                'GC_Content': fwd['GC_Content']
                            })

        # Save the valid outer primer combinations if any are found
        if outer_primer_combinations:
            outer_primer_combinations_df = pd.DataFrame(outer_primer_combinations)
            outer_primer_combinations_df.drop_duplicates(subset=['Combination_Name'], inplace=True)

            outer_output_file = 'outer_semi_nested_primer_combinations.tsv'
            with open(outer_output_file, 'w') as f:
                # Add header
                f.write("Forward_Sequence\tReverse_Sequence\tCombination_Name\n")
                for idx, row in outer_primer_combinations_df.iterrows():
                    f.write(f"{row['Forward_Sequence']}\t{row['Reverse_Sequence']}\t{row['Combination_Name']}\n")

            logging.info(f"Formatted outer primer combinations saved to {outer_output_file}")
        else:
            logging.info("No valid outer primer combinations found.")

        logging.info("Searching for valid semi-nested primer combinations")
        for outer_comb in tqdm(outer_primer_combinations, desc="Processing Semi-Nested Primers"):
            fwd_outer_start = outer_comb['Forward_Start']
            rev_outer_end = outer_comb['Reverse_End']

            inner_fwd_primers = fwd_primers[(fwd_primers['Start'] > fwd_outer_start + 50) & 
                                            (fwd_primers['Start'] < rev_outer_end - 50)]
            inner_rev_primers = rev_primers[(rev_primers['Start'] > fwd_outer_start + 50) & 
                                            (rev_primers['Start'] < rev_outer_end - 50)]

            # Use the outer forward primer with an inner reverse primer
            for _, rev in inner_rev_primers.iterrows():
                if not check_primer_conditions(rev):
                    continue

                combination_name = f"{outer_comb['Forward_Primer']}_{rev['Primer']}"
                tm_max_diff = abs(rev['Tm_max'] - outer_comb['Tm_max'])
                tm_min_diff = abs(rev['Tm_min'] - outer_comb['Tm_min'])
                gc_content_diff = abs(rev['GC_Content'] - outer_comb['GC_Content'])

                logging.debug(f"Evaluating semi-nested combination: {combination_name} | Tm Max Diff: {tm_max_diff}, Tm Min Diff: {tm_min_diff}, GC Content Diff: {gc_content_diff}")

                if (tm_max_diff <= 3) and (tm_min_diff <= 3) and (gc_content_diff <= 10):
                    semi_nested_primer_combinations.append({
                        'Forward_Primer': outer_comb['Forward_Sequence'],
                        'Forward_Sequence': outer_comb['Forward_Sequence'],
                        'Reverse_Primer': rev['Primer'],
                        'Reverse_Sequence': rev['Sequence'],
                        'Combination_Name': combination_name.replace(' ', '_'),
                        'Outer_Combination': outer_comb['Combination_Name']
                    })

            # Use the outer reverse primer with an inner forward primer
            for _, fwd in inner_fwd_primers.iterrows():
                if not check_primer_conditions(fwd):
                    continue

                combination_name = f"{fwd['Primer']}_{outer_comb['Reverse_Primer']}"
                tm_max_diff = abs(fwd['Tm_max'] - outer_comb['Tm_max'])
                tm_min_diff = abs(fwd['Tm_min'] - outer_comb['Tm_min'])
                gc_content_diff = abs(fwd['GC_Content'] - outer_comb['GC_Content'])

                logging.debug(f"Evaluating semi-nested combination: {combination_name} | Tm Max Diff: {tm_max_diff}, Tm Min Diff: {tm_min_diff}, GC Content Diff: {gc_content_diff}")

                if (tm_max_diff <= 3) and (tm_min_diff <= 3) and (gc_content_diff <= 10):
                    semi_nested_primer_combinations.append({
                        'Forward_Primer': fwd['Primer'],
                        'Forward_Sequence': fwd['Sequence'],
                        'Reverse_Primer': outer_comb['Reverse_Primer'],
                        'Reverse_Sequence': outer_comb['Reverse_Sequence'],
                        'Combination_Name': combination_name.replace(' ', '_'),
                        'Outer_Combination': outer_comb['Combination_Name']
                    })

        # Save the valid semi-nested primer combinations if any are found
        if semi_nested_primer_combinations:
            semi_nested_primer_combinations_df = pd.DataFrame(semi_nested_primer_combinations)
            semi_nested_primer_combinations_df.drop_duplicates(subset=['Combination_Name'], inplace=True)
            
            semi_nested_output_file = 'inner_semi_nested_primer_combinations.tsv'
            with open(semi_nested_output_file, 'w') as f:
                # Add header
                f.write("Forward_Sequence\tReverse_Sequence\tCombination_Name\tOuter_Combination\n")
                for idx, row in semi_nested_primer_combinations_df.iterrows():
                    f.write(f"{row['Forward_Sequence']}\t{row['Reverse_Sequence']}\t{row['Combination_Name']}\t{row['Outer_Combination']}\n")
            
            logging.info(f"Formatted semi-nested primer combinations saved to {semi_nested_output_file}")
        else:
            logging.info("No valid semi-nested primer combinations found.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python generate_semi_nested_primer_combinations.py mapping_positions.tsv")
        sys.exit(1)
    else:
        setup_logging()
        input_file = sys.argv[1]
        regions_of_interest = {
            'NS1': (1000, 2000),
            'NS3': (3000, 4000),
            'NS5': (5000, 7000)
        }
        generate_semi_nested_primer_combinations(input_file, regions_of_interest)
