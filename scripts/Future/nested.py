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

def generate_nested_primer_combinations(input_file, regions_of_interest):
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
        inner_primer_combinations = []

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
                        if fwd['Start'] >= start and rev['Start'] <= end:
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

            outer_output_file = 'outer_nested_primer_combinations.tsv'
            with open(outer_output_file, 'w') as f:
                # Add header
                f.write("Forward_Sequence\tReverse_Sequence\tCombination_Name\n")
                for idx, row in outer_primer_combinations_df.iterrows():
                    f.write(f"{row['Forward_Sequence']}\t{row['Reverse_Sequence']}\t{row['Combination_Name']}\n")

            logging.info(f"Formatted outer primer combinations saved to {outer_output_file}")
        else:
            logging.info("No valid outer primer combinations found.")

        logging.info("Searching for valid inner primer combinations")
        for outer_comb in tqdm(outer_primer_combinations, desc="Processing Inner Primers"):
            fwd_outer_start = outer_comb['Forward_Start']
            rev_outer_end = outer_comb['Reverse_End']

            inner_fwd_primers = fwd_primers[(fwd_primers['Start'] > fwd_outer_start + 50) & 
                                            (fwd_primers['Start'] < rev_outer_end - 50)]
            inner_rev_primers = rev_primers[(rev_primers['Start'] > fwd_outer_start + 50) & 
                                            (rev_primers['Start'] < rev_outer_end - 50)]

            for _, fwd in inner_fwd_primers.iterrows():
                if not check_primer_conditions(fwd):
                    continue

                matching_revs = inner_rev_primers[(inner_rev_primers['Reference'] == fwd['Reference']) & 
                                                  ((inner_rev_primers['Genotype'] == fwd['Genotype']) | 
                                                   (inner_rev_primers['Genotype'] == 'ALL')) &
                                                  (inner_rev_primers['Start'] > fwd['Start'] + 100) & 
                                                  (inner_rev_primers['Start'] < fwd['Start'] + 500)]
                matching_revs['Amplicon_Length'] = matching_revs['Start'] - fwd['Start']

                for _, rev in matching_revs.iterrows():
                    if not check_primer_conditions(rev):
                        continue

                    combination_name = f"{fwd['Primer']}_{rev['Primer']}"
                    tm_max_diff = abs(fwd['Tm_max'] - rev['Tm_max'])
                    tm_min_diff = abs(fwd['Tm_min'] - rev['Tm_min'])
                    gc_content_diff = abs(fwd['GC_Content'] - rev['GC_Content'])

                    if (tm_max_diff <= 3) and (tm_min_diff <= 3) and (gc_content_diff <= 10):
                        combination_name = f"{fwd['Primer']}_{rev['Primer']}"
                        inner_primer_combinations.append({
                            'Forward_Primer': fwd['Primer'],
                            'Forward_Sequence': fwd['Sequence'],
                            'Reverse_Primer': rev['Primer'],
                            'Reverse_Sequence': rev['Sequence'],
                            'Combination_Name': combination_name.replace(' ', '_'),
                            'Outer_Combination': outer_comb['Combination_Name'],
                            'Amplicon_Length': rev['Amplicon_Length'],
                            'Tm_max_diff': tm_max_diff,
                            'Tm_min_diff': tm_min_diff,
                            'GC_Content_diff': gc_content_diff
                        })
                        logging.debug(f"Valid inner combination found: {combination_name} with Amplicon Length: {rev['Amplicon_Length']}, Tm Max Diff: {tm_max_diff}, Tm Min Diff: {tm_min_diff}, GC Content Diff: {gc_content_diff}")

        # Save the valid inner primer combinations if any are found
        if inner_primer_combinations:
            inner_primer_combinations_df = pd.DataFrame(inner_primer_combinations)
            inner_primer_combinations_df.drop_duplicates(subset=['Combination_Name'], inplace=True)
            
            inner_output_file = 'inner_nested_primer_combinations.tsv'
            with open(inner_output_file, 'w') as f:
                # Add header
                f.write("Forward_Sequence\tReverse_Sequence\tCombination_Name\tOuter_Combination\n")
                for idx, row in inner_primer_combinations_df.iterrows():
                    f.write(f"{row['Forward_Sequence']}\t{row['Reverse_Sequence']}\t{row['Combination_Name']}\t{row['Outer_Combination']}\n")

            logging.info(f"Formatted inner primer combinations saved to {inner_output_file}")
        else:
            logging.info("No valid inner primer combinations found.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python generate_nested_primer_combinations.py mapping_positions.tsv")
        sys.exit(1)
    else:
        setup_logging()
        input_file = sys.argv[1]
        regions_of_interest = {
            'NS1': (1000, 2000),
            'NS3': (3000, 4000),
            'NS5': (5000, 7000)
        }
        generate_nested_primer_combinations(input_file, regions_of_interest)
