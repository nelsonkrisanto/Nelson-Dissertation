import pandas as pd
import sys

def generate_nested_primer_combinations(input_file):
    try:
        print(f"Loading mapping data from {input_file}")
        df_mapping = pd.read_csv(input_file, delimiter='\t')

        print("Loading primer metadata from primer_metadata.tsv")
        df_metadata = pd.read_csv('primer_metadata.tsv', delimiter='\t')

        print("Merging dataframes based on 'Primer' column")
        merged_df = pd.merge(df_mapping, df_metadata, on='Primer', how='inner')
        print("Merged dataframes successfully.")

        primer_combinations_first_round = []

        print("Searching for first round valid primer combinations")
        for _, row1 in merged_df.iterrows():
            if row1['Primer'].endswith('_F'):
                for _, row2 in merged_df.iterrows():
                    if row2['Primer'].endswith('_R') and row1['Reference'] == row2['Reference']:
                        diff = abs(row1['Start'] - row2['Start'])
                        if diff > 150:
                            amplicon_length = row2['Start'] - row1['Start']
                            if 100 < amplicon_length < 1000:  # Filter amplicon length for the first round
                                combination_name = f"{row1['Primer']}_{row2['Primer']}"
                                tm_max_diff = round(abs(row1['Tm_max'] - row2['Tm_max']))
                                tm_min_diff = round(abs(row1['Tm_min'] - row2['Tm_min']))

                                if tm_max_diff <= 7 and tm_min_diff <= 7:
                                    primer_combinations_first_round.append({
                                        'Forward_Primer': row1['Primer'],
                                        'Reverse_Primer': row2['Primer'],
                                        'Combination_Name': combination_name,
                                        'Reference': row1['Reference'],
                                        'Amplicon_Length_First_Round': amplicon_length
                                    })
                                    print(f"Found a valid first-round primer combination: {combination_name}")

        # Now, find valid second-round primer combinations within the first-round amplicons
        primer_combinations_second_round = []

        print("Searching for second round valid primer combinations within first-round amplicons")
        for combo in primer_combinations_first_round:
            forward_primers = merged_df[(merged_df['Primer'].endswith('_F')) & (merged_df['Reference'] == combo['Reference'])]
            reverse_primers = merged_df[(merged_df['Primer'].endswith('_R')) & (merged_df['Reference'] == combo['Reference'])]
            
            for _, fwd in forward_primers.iterrows():
                for _, rev in reverse_primers.iterrows():
                    if fwd['Start'] > combo['Forward_Primer_Start'] and rev['Start'] < combo['Reverse_Primer_Start']:
                        amplicon_length_second = rev['Start'] - fwd['Start']
                        if 80 < amplicon_length_second < 500:  # Nested amplicon length range
                            second_round_name = f"{fwd['Primer']}-{rev['Primer']} within {combo['Combination_Name']}"
                            primer_combinations_second_round.append({
                                'First_Round_Combination': combo['Combination_Name'],
                                'Second_Round_Forward_Primer': fwd['Primer'],
                                'Second_Round_Reverse_Primer': rev['Primer'],
                                'Second_Round_Combination_Name': second_round_name,
                                'Amplicon_Length_Second_Round': amplicon_length_second
                            })
                            print(f"Found a valid second-round primer combination: {second_round_name}")

        # Save the results
        if primer_combinations_second_round:
            try:
                primer_combinations_df = pd.DataFrame(primer_combinations_second_round)
                primer_combinations_df.drop_duplicates(subset='Second_Round_Combination_Name', inplace=True)

                output_file = 'nested_primer_combinations.tsv'
                primer_combinations_df.to_csv(output_file, sep='\t', index=False)
                print(f"Nested primer combinations saved to {output_file}")

            except Exception as e:
                print(f"An error occurred while saving nested primer combinations: {e}")
        else:
            print("No valid nested primer combinations found.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python nested_combine.py mapping_positions.tsv")
    else:
        input_file = sys.argv[1]
        generate_nested_primer_combinations(input_file)
