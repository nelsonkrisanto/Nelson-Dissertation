import pandas as pd
import sys

def calculate_average_tm(min_tm, max_tm):
    avg_tm = (min_tm + max_tm) / 2
    return avg_tm

def check_Tm_condition(tm, min_tm, max_tm):
    avg_tm = calculate_average_tm(min_tm, max_tm)
    return (tm >= min_tm and tm <= max_tm) or abs(tm - min_tm) <= 10 or abs(max_tm - tm) <= 10

def main(mapping_positions_file, primer_metadata_file, primer_combinations_file):
    # Read data from files
    mapping_positions = pd.read_csv(mapping_positions_file, sep="\t")
    primer_metadata = pd.read_csv(primer_metadata_file, sep="\t")
    primer_combinations = pd.read_csv(primer_combinations_file, sep="\t")

    # Merge mapping_positions and primer_metadata
    mapping_positions = mapping_positions.merge(primer_metadata, on="Primer")

    result_df = pd.DataFrame(columns=["Combination_Name", "Valid_Primer_Count"])
    selected_primer_info = []
    unique_pcr_counter = 1

    for i, row in primer_combinations.iterrows():
        combination_name = row["Combination_Name"]
        forward_start = row["Forward_Start"]
        reverse_start = row["Reverse_Start"]
        genogroup = row["Genogroup"]
        ref = row["Reference"]

    # Find the corresponding forward and reverse primers in mapping_positions
    forward_primer = row["Forward_Primer"]
    reverse_primer = row["Reverse_Primer"]
    
    # Extract Tm values for the forward primer from mapping_positions
    tm_min = mapping_positions.loc[mapping_positions["Primer"] == forward_primer, "Tm_min"].values[0]
    tm_max = mapping_positions.loc[mapping_positions["Primer"] == forward_primer, "Tm_max"].values[0]
    
    print(f"Tm values for {forward_primer}: min = {tm_min}, max = {tm_max}")
    
    tm = calculate_average_tm(tm_min, tm_max)
    
    # Filter reverse primers that meet the criteria
    valid_reverse_primers = reverse_primers[
        (reverse_primers["Reference"] == ref) &
        (reverse_primers["Start"] > (forward_start + 100)) &
        (reverse_primers["Start"] < reverse_start)
    ]
    
    print(f"Valid reverse primers before temperature check: {len(valid_reverse_primers)}")
    
    # Filter reverse primers based on temperature conditions
    valid_reverse_primers = valid_reverse_primers[
        valid_reverse_primers.apply(lambda row: check_Tm_condition(tm, tm_min, tm_max), axis=1)
    ]
    
    print(f"Valid reverse primers after temperature check: {len(valid_reverse_primers)}")
    
    # Filter valid_reverse_primers using the valid primer names and other conditions
    valid_reverse_primers = valid_reverse_primers[
        (~valid_reverse_primers["Combination_Name"].str.contains(forward_primer)) &
        (valid_reverse_primers["Genogroup"] == genogroup) &
        (valid_reverse_primers["Reference"] == ref)
    ]
    
    print(f"Filtered reverse primers: {len(valid_reverse_primers)}")
    
    # Store the result in the result_df
    result_df = result_df.append({"Combination_Name": combination_name, "Valid_Primer_Count": len(valid_reverse_primers)}, ignore_index=True)
    
    if len(valid_reverse_primers) > 0:
        reverse_end = valid_reverse_primers.loc[valid_reverse_primers.index[0], "End"]
        
        # Calculate amplicon lengths
        first_round_amplicon = row["Amplicon_Length"]
        second_round_amplicon = reverse_end - forward_start
        
        for j in range(len(first_round_amplicon)):
            if first_round_amplicon[j] > second_round_amplicon[j]:
                print(f"Storing selected information for combination: {combination_name}")
                
                selected_primer_info.append({
                    "first_round_pcr": combination_name,
                    "second_round_pcr": f"{forward_primer}-{valid_reverse_primers.loc[valid_reverse_primers.index[0], 'Primer']}",
                    "nested_pcr": f"{combination_name}_{forward_primer}_{valid_reverse_primers.loc[valid_reverse_primers.index[0], 'Primer']}",
                    "Forward_Start": forward_start,
                    "Reverse_Start": valid_reverse_primers.loc[valid_reverse_primers.index[0], "Start"],
                    "first_round_amplicon": first_round_amplicon[j],
                    "second_round_amplicon": second_round_amplicon[j],
                    "Genogroup": genogroup,
                    "Reference": ref,
                    "New_Primer": valid_reverse_primers.loc[valid_reverse_primers.index[0], "Primer"],
                    "nested_pcr_type": "Semi_Nested_PCR"
                })
                
                unique_pcr_counter += 1
            else:
                print(f"Amplicon length condition not satisfied for combination: {combination_name}")

# Create a DataFrame from selected_primer_info list
    result_df_combined = pd.DataFrame(selected_primer_info)

    # Drop duplicates from the combined result DataFrame
    result_df_combined = result_df_combined.drop_duplicates()

    # Filter result_df_combined to only include rows with second_round_amplicon > 100
    result_df_filtered = result_df_combined[result_df_combined["second_round_amplicon"] > 100]

    # Export the DataFrames to TSV files
    result_df_filtered.to_csv("semi_nested_reverse.tsv", sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python program_name.py mapping_positions.tsv primer_metadata.tsv primer_combinations.tsv")
    else:
        mapping_positions_file = sys.argv[1]
        primer_metadata_file = sys.argv[2]
        primer_combinations_file = sys.argv[3]
        main(mapping_positions_file, primer_metadata_file, primer_combinations_file)
