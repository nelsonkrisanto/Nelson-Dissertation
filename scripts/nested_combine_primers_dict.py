import pandas as pd
import sys
import logging
from tqdm import tqdm

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def load_results(outer_file, inner_file):
    logging.info(f"Loading outer results from {outer_file}")
    outer_df = pd.read_csv(outer_file, delimiter='\t', header=None, names=['Forward_Primer', 'Reverse_Primer', 'Combination_Name'])
    logging.info(f"Loading inner results from {inner_file}")
    inner_df = pd.read_csv(inner_file, delimiter='\t', header=None, names=['Forward_Primer', 'Reverse_Primer', 'Combination_Name'])
    return outer_df, inner_df

def map_inner_to_outer(outer_df, inner_df):
    combined_results = []
    logging.info("Mapping inner to outer results...")

    # Create a dictionary for quick lookup
    outer_dict = {outer_row['Combination_Name']: outer_row for _, outer_row in outer_df.iterrows()}

    for inner_idx, inner_row in tqdm(inner_df.iterrows(), total=inner_df.shape[0], desc="Processing inner primers"):
        if isinstance(inner_row['Combination_Name'], str):
            # Get the outer combination name part to match
            outer_name_part = inner_row['Combination_Name'].split('_')[0]
            if outer_name_part in outer_dict:
                outer_row = outer_dict[outer_name_part]
                combined_results.append({
                    'Outer_Forward_Primer': outer_row['Forward_Primer'],
                    'Outer_Reverse_Primer': outer_row['Reverse_Primer'],
                    'Outer_Combination_Name': outer_row['Combination_Name'],
                    'Inner_Forward_Primer': inner_row['Forward_Primer'],
                    'Inner_Reverse_Primer': inner_row['Reverse_Primer'],
                    'Inner_Combination_Name': inner_row['Combination_Name']
                })

    combined_df = pd.DataFrame(combined_results)
    logging.info(f"Total combined results: {len(combined_df)}")
    return combined_df

def save_combined_results(combined_df, output_file):
    combined_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Combined results saved to {output_file}")

if __name__ == "__main__":
    setup_logging()
    
    if len(sys.argv) != 4:
        print("Usage: python combine_nested_insilico_results.py outer_results_file inner_results_file combined_output_file")
        sys.exit(1)

    outer_results_file = sys.argv[1]
    inner_results_file = sys.argv[2]
    combined_output_file = sys.argv[3]

    logging.info("Starting the combination process...")
    outer_df, inner_df = load_results(outer_results_file, inner_results_file)
    combined_df = map_inner_to_outer(outer_df, inner_df)
    save_combined_results(combined_df, combined_output_file)
    logging.info("Combination process completed successfully.")
