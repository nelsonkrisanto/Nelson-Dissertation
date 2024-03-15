import sys
import pandas as pd

def calculate_average_tm(min_tm, max_tm):
    avg_tm = (min_tm + max_tm) / 2
    return avg_tm

def check_Tm_condition(tm, min_tm, max_tm):
    avg_tm = calculate_average_tm(min_tm, max_tm)
    return (tm >= min_tm and tm <= max_tm) or abs(tm - min_tm) <= 10 or abs(max_tm - tm) <= 10

def main(mapping_positions_path, primer_metadata_path, primer_combinations_path):
    # Read mapping_positions and primer_combinations data
    mapping_positions = pd.read_table(mapping_positions_path, sep="\t")
    primer_metadata = pd.read_table(primer_metadata_path, sep="\t")
    primer_combinations = pd.read_table(primer_combinations_path, sep="\t")

    # Merge mapping_positions with primer_metadata
    mapping_positions = mapping_positions.merge(primer_metadata, on="Primer")

    # Filter forward primers from mapping_positions
    forward_primers = mapping_positions[mapping_positions["Primer"].str.endswith("F")]

    # Initialize a counter for unique PCR names
    unique_pcr_counter = 1

    selected_primer_info = []

    for i in range(len(primer_combinations)):
        combination_name = primer_combinations.loc[i, "Combination_Name"]
        forward_start = primer_combinations.loc[i, "Forward_Start"]
        reverse_start = primer_combinations.loc[i, "Reverse_Start"]
        genogroup = primer_combinations.loc[i, "Genogroup"]
        ref = primer_combinations.loc[i, "Reference"]
        
        forward_primer = primer_combinations.loc[i, "Forward_Primer"]
        reverse_primer = primer_combinations.loc[i, "Reverse_Primer"]
        
        tm_min = mapping_positions.loc[mapping_positions["Primer"] == forward_primer, "Tm_min"].values[0]
        tm_max = mapping_positions.loc[mapping_positions["Primer"] == forward_primer, "Tm_max"].values[0]
        
        tm = calculate_average_tm(tm_min, tm_max)
        
        valid_forward_primers = forward_primers[
            (forward_primers["Reference"] == ref) &
            (forward_primers["Start"] > (reverse_start - 100)) &
            (forward_primers["Start"] > forward_start)
        ]
        
        valid_forward_primers = valid_forward_primers[
            valid_forward_primers.apply(lambda row: check_Tm_condition(tm, tm_min, tm_max), axis=1)
        ]
        
        valid_forward_primers = valid_forward_primers[
            (~valid_forward_primers["Combination_Name"].str.contains(forward_primer)) &
            (valid_forward_primers["Genogroup"] == genogroup) &
            (valid_forward_primers["Reference"] == ref)
        ]
        
        if len(valid_forward_primers) > 0:
            reverse_end = mapping_positions[mapping_positions["Primer"] == valid_forward_primers.iloc[0]["Primer"]]["End"].values[0]
            
            for first_round_amplicon, second_round_amplicon in zip(primer_combinations.loc[i, "Amplicon_Length"], reverse_end - forward_start):
                if first_round_amplicon > second_round_amplicon:
                    selected_primer_info.append({
                        "first_round_pcr": combination_name,
                        "second_round_pcr": f"{valid_forward_primers.iloc[0]['Primer']}-{reverse_primer}",
                        "nested_pcr": f"{combination_name}_{valid_forward_primers.iloc[0]['Primer']}_{reverse_primer}",
                        "Forward_Start": valid_forward_primers.iloc[0]["Start"],
                        "Reverse_Start": reverse_start,
                        "first_round_amplicon": first_round_amplicon,
                        "second_round_amplicon": second_round_amplicon,
                        "Genogroup": genogroup,
                        "Reference": ref,
                        "New_Primer": valid_forward_primers.iloc[0]["Primer"],
                        "nested_pcr_type": "Semi_Nested_PCR"
                    })
                    unique_pcr_counter += 1
    
    # Create a DataFrame from the selected_primer_info list
    result_df_combined = pd.DataFrame(selected_primer_info)

    # Drop duplicate rows from result_df_combined
    result_df_combined = result_df_combined.drop_duplicates()

    # Filter rows with second_round_amplicon > 100
    result_df_filtered = result_df_combined[result_df_combined["second_round_amplicon"] > 100]

    # Save the filtered result to a TSV file named semi_nested_forward.tsv
    result_df_filtered.to_csv("semi_nested_forward.tsv", sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script_name.py mapping_positions.tsv primer_metadata.tsv primer_combinations.tsv")
    else:
        mapping_positions_path = sys.argv[1]
        primer_metadata_path = sys.argv[2]
        primer_combinations_path = sys.argv[3]
        main(mapping_positions_path, primer_metadata_path, primer_combinations_path)
