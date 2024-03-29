import pandas as pd
import sys
from tqdm import tqdm

def generate_primer_combinations(input_file):
    try:
        print(f"Loading mapping data from {input_file}")
        df_mapping = pd.read_csv(input_file, delimiter='\t')

        print("Loading primer metadata from primer_metadata.tsv")
        df_metadata = pd.read_csv('primer_metadata.tsv', delimiter='\t')

        print("Merging dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        print("Merged dataframes successfully.")

        # Separate forward and reverse primers into two dataframes
        fwd_primers = merged_df[merged_df['Primer'].str.endswith('_F')]
        rev_primers = merged_df[merged_df['Primer'].str.endswith('_R')]

        primer_combinations = []

        print("Searching for valid primer combinations")
        for _, fwd in tqdm(fwd_primers.iterrows(), total=fwd_primers.shape[0], desc="Processing Primers"):
            # Filter reverse primers from the same reference and where start position difference is > 150
            matching_revs = rev_primers[(rev_primers['Reference'] == fwd['Reference']) & 
                                        (rev_primers['Start'] > fwd['Start'] + 150) & 
                                        (rev_primers['Start'] < fwd['Start'] + 1000)]
            
            for _, rev in matching_revs.iterrows():
                combination_name = f"{fwd['Primer']}_{rev['Primer']}"
                amplicon_length = rev['Start'] - fwd['Start']
                tm_max_diff = round(abs(fwd['Tm_max'] - rev['Tm_max']))
                tm_min_diff = round(abs(fwd['Tm_min'] - rev['Tm_min']))

                if tm_max_diff <= 7 and tm_min_diff <= 7:
                    primer_combinations.append({
                        'Forward_Primer': fwd['Primer'],
                        'Reverse_Primer': rev['Primer'],
                        'Combination_Name': combination_name,
                        'Reference': fwd['Reference'],
                        'Amplicon_Length': amplicon_length
                    })

        if primer_combinations:
            primer_combinations_df = pd.DataFrame(primer_combinations)
            primer_combinations_df.drop_duplicates(subset='Combination_Name', inplace=True)

            output_file = 'primer_combinations.tsv'
            primer_combinations_df.to_csv(output_file, sep='\t', index=False)
            print(f"Primer combinations saved to {output_file}")
        else:
            print("No valid primer combinations found.")

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python combine.py mapping_positions.tsv")
        sys.exit(1)
    else:
        input_file = sys.argv[1]
        generate_primer_combinations(input_file)
