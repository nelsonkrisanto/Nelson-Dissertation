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

def load_metadata(mapping_file, metadata_file):
    logging.info(f"Loading mapping data from {mapping_file}")
    df_mapping = pd.read_csv(mapping_file, delimiter='\t')
    logging.info(f"Mapping data loaded: {df_mapping.shape[0]} rows")
    
    logging.info(f"Loading primer metadata from {metadata_file}")
    df_metadata = pd.read_csv(metadata_file, delimiter='\t')
    logging.info(f"Primer metadata loaded: {df_metadata.shape[0]} rows")
    
    return df_mapping, df_metadata

def add_genotype_reference_info(df, mapping_df, metadata_df):
    merged_df = pd.merge(df, mapping_df[['Primer', 'Genotype', 'Reference']], left_on='Forward_Primer', right_on='Primer', how='left')
    merged_df = pd.merge(merged_df, metadata_df[['Primer', 'Genotype']], left_on='Forward_Primer', right_on='Primer', how='left', suffixes=('_map', '_meta'))
    merged_df['Genotype'] = merged_df['Genotype_map'].combine_first(merged_df['Genotype_meta'])
    merged_df.drop(columns=['Primer_map', 'Primer_meta'], inplace=True)
    
    missing_genotype = merged_df['Genotype'].isnull().sum()
    missing_reference = merged_df['Reference'].isnull().sum()
    
    logging.debug(f"Missing Genotype after merge: {missing_genotype}")
    logging.debug(f"Missing Reference after merge: {missing_reference}")
    
    return merged_df

def map_inner_to_outer(outer_df, inner_df):
    combined_results = []
    logging.info("Mapping inner to outer results...")

    outer_dict = {(outer_row['Combination_Name'], outer_row['Genotype'], outer_row['Reference']): outer_row for _, outer_row in outer_df.iterrows()}

    for inner_idx, inner_row in tqdm(inner_df.iterrows(), total=inner_df.shape[0], desc="Processing inner primers"):
        if isinstance(inner_row['Combination_Name'], str):
            outer_name_part = inner_row['Combination_Name'].split('_')[0]
            genotype = inner_row['Genotype']
            reference = inner_row['Reference']
            key = (outer_name_part, genotype, reference)
            
            if key in outer_dict:
                outer_row = outer_dict[key]
                combined_results.append({
                    'Outer_Forward_Primer': outer_row['Forward_Primer'],
                    'Outer_Reverse_Primer': outer_row['Reverse_Primer'],
                    'Outer_Combination_Name': outer_row['Combination_Name'],
                    'Outer_Genotype': outer_row['Genotype'],
                    'Outer_Reference': outer_row['Reference'],
                    'Inner_Forward_Primer': inner_row['Forward_Primer'],
                    'Inner_Reverse_Primer': inner_row['Reverse_Primer'],
                    'Inner_Combination_Name': inner_row['Combination_Name'],
                    'Inner_Genotype': inner_row['Genotype'],
                    'Inner_Reference': inner_row['Reference']
                })
            else:
                logging.warning(f"No match for inner primer combination: {inner_row['Combination_Name']} with genotype: {genotype} and reference: {reference}")

    combined_df = pd.DataFrame(combined_results)
    logging.info(f"Total combined results: {len(combined_df)}")
    return combined_df

def save_combined_results(combined_df, output_file):
    combined_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Combined results saved to {output_file}")

if __name__ == "__main__":
    setup_logging()
    
    if len(sys.argv) != 6:
        print("Usage: python combine_nested_insilico_results.py outer_results_file inner_results_file combined_output_file mapping_file metadata_file")
        sys.exit(1)

    outer_results_file = sys.argv[1]
    inner_results_file = sys.argv[2]
    combined_output_file = sys.argv[3]
    mapping_file = sys.argv[4]
    metadata_file = sys.argv[5]

    logging.info("Starting the combination process...")
    outer_df, inner_df = load_results(outer_results_file, inner_results_file)
    mapping_df, metadata_df = load_metadata(mapping_file, metadata_file)
    
    outer_df = add_genotype_reference_info(outer_df, mapping_df, metadata_df)
    inner_df = add_genotype_reference_info(inner_df, mapping_df, metadata_df)
    
    combined_df = map_inner_to_outer(outer_df, inner_df)
    save_combined_results(combined_df, combined_output_file)
    logging.info("Combination process completed successfully.")
