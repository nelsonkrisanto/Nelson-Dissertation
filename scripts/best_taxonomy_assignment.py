import sys
import pandas as pd

def main(taxonomy_file_path, output_file_path):
    # Load the combined taxonomy assignment results
    taxonomy_data = pd.read_csv(taxonomy_file_path, sep='\t', header=None)

    # Split the data into columns
    taxonomy_data[['Sequence_ID', 'Taxonomy', 'Metadata']] = taxonomy_data[0].str.split(pat=' ', n=2, expand=True)

    # Split the Metadata column to extract percent identity
    taxonomy_data[['Percent_ID', 'Length_Hit', 'Length_Query']] = taxonomy_data['Metadata'].str.split(pat=':', expand=True)

    # Convert Percent_ID to numeric for sorting
    taxonomy_data['Percent_ID'] = pd.to_numeric(taxonomy_data['Percent_ID'])

    # Filter out entries that do not meet the minimum cutoff for phylum (85%)
    taxonomy_data = taxonomy_data[taxonomy_data['Percent_ID'] >= 65]

    # Select the taxonomy with the highest percent identity for each sequence
    best_taxonomy = taxonomy_data.loc[taxonomy_data.groupby('Sequence_ID')['Percent_ID'].idxmax()]

    # Drop unnecessary columns
    best_taxonomy = best_taxonomy.drop(columns=[0, 'Metadata', 'Length_Hit', 'Length_Query'])

    # Save the best taxonomy assignment to a new file
    best_taxonomy.to_csv(output_file_path, sep='\t', index=False, header=False)

    print(f"Best taxonomy assignment per sequence has been saved to {output_file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python best_taxonomy_assignment.py <taxonomy_assignment_file> <output_file>")
        sys.exit(1)

    taxonomy_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    main(taxonomy_file_path, output_file_path)
