import pandas as pd

def calculate_average_tm(min_tm, max_tm):
    avg_tm = (min_tm + max_tm) / 2
    return avg_tm

def check_tm_condition(tm, min_tm, max_tm):
    avg_tm = calculate_average_tm(min_tm, max_tm)
    return (tm >= min_tm and tm <= max_tm) or abs(tm - min_tm) <= 10 or abs(max_tm - tm) <= 10

def main():
    mapping_positions_file = input("Enter the path to mapping_positions.tsv: ")
    primer_metadata_file = input("Enter the path to primer_metadata.tsv: ")
    primer_combinations_file = input("Enter the path to primer_combinations.tsv: ")

    mapping_positions = pd.read_table(mapping_positions_file, sep="\t")
    primer_metadata = pd.read_table(primer_metadata_file, sep="\t")
    primer_combinations = pd.read_table(primer_combinations_file, sep="\t")

    mapping_positions = mapping_positions.merge(primer_metadata, on='Primer', how='left')
    reverse_primers = mapping_positions[mapping_positions['Primer'].str.endswith('R')]

    result_df = pd.DataFrame(columns=['Combination_Name', 'Valid_Primer_Count'])
    selected_primer_info = []
    unique_pcr_counter = 1

    for i in range(len(primer_combinations)):
        combination_name = primer_combinations['Combination_Name'][i]
        forward_start = primer_combinations['Forward_Start'][i]
        reverse_start = primer_combinations['Reverse_Start'][i]
        genogroup = primer_combinations['Genogroup'][i]
        ref = primer_combinations['Reference'][i]
        
        forward_primer = primer_combinations['Forward_Primer'][i]
        reverse_primer = primer_combinations['Reverse_Primer'][i]
        
        tm_min = mapping_positions.loc[mapping_positions['Primer'] == forward_primer, 'Tm_min'].values[0]
        tm_max = mapping_positions.loc[mapping_positions['Primer'] == forward_primer, 'Tm_max'].values[0]
        
        tm = calculate_average_tm(tm_min, tm_max)
        
        valid_reverse_primers = reverse_primers[
            (reverse_primers['Reference'] == ref) &
            (reverse_primers['Start'] > (forward_start + 100)) &
            (reverse_primers['Start'] < reverse_start)
        ]
        
        valid_reverse_primers = valid_reverse_primers[
            valid_reverse_primers.apply(lambda row: check_tm_condition(tm, tm_min, tm_max), axis=1)
        ]
        
        valid_reverse_primers = valid_reverse_primers[
            (~valid_reverse_primers['Combination_Name'].str.contains(forward_primer)) &
            (valid_reverse_primers['Genogroup'] == genogroup) &
            (valid_reverse_primers['Reference'] == ref)
        ]
        
        result_df = result_df.append({
            'Combination_Name': combination_name,
            'Valid_Primer_Count': len(valid_reverse_primers)
        }, ignore_index=True)
        
        if len(valid_reverse_primers) > 0:
            reverse_end = reverse_primers[reverse_primers['Primer'] == valid_reverse_primers['Primer'].iloc[0]]['End'].values[0]
            
            first_round_amplicon = primer_combinations['Amplicon_Length'][i]
            second_round_amplicon = reverse_end - forward_start
            
            for j in range(len(first_round_amplicon)):
                if first_round_amplicon[j] > second_round_amplicon[j]:
                    print(f"Storing selected information for combination: {combination_name}")
                    selected_primer_info.append({
                        'first_round_pcr': combination_name,
                        'second_round_pcr': f"{primer_combinations['Forward_Primer'][i]}-{valid_reverse_primers['Primer'].iloc[0]}",
                        'nested_pcr': f"{combination_name}_{primer_combinations['Forward_Primer'][i]}_{valid_reverse_primers['Primer'].iloc[0]}",
                        'Forward_Start': forward_start,
                        'Reverse_Start': valid_reverse_primers['Start'].iloc[0],
                        'first_round_amplicon': first_round_amplicon[j],
                        'second_round_amplicon': second_round_amplicon[j],
                        'Genogroup': genogroup,
                        'Reference': ref,
                        'New_Primer': valid_reverse_primers['Primer'].iloc[0],
                        'nested_pcr_type': "Semi_Nested_PCR"
                    })
                    unique_pcr_counter += 1
                else:
                    print(f"Amplicon length condition not satisfied for combination: {combination_name}")

    result_df_combined = pd.DataFrame(selected_primer_info)
    result_df_combined = result_df_combined.drop_duplicates()
    result_df_filtered = result_df_combined[result_df_combined['second_round_amplicon'] > 100]

    result_df_filtered.to_csv("nested_primers.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
