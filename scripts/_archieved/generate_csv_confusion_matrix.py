import pandas as pd
import sys

def generate_final_csv(amplicons_genotyping_file, whole_genome_genotyping_file, primer_file, output_csv):
    # Load the taxonomy assignment results for amplicons and whole genome
    df_amplicons = pd.read_csv(amplicons_genotyping_file, sep='\t', header=None, names=[
        'Amplicons', 'Taxonomy', 'Score'])
    df_whole_genome = pd.read_csv(whole_genome_genotyping_file, sep='\t', header=None, names=[
        'Accession_Number', 'Taxonomy', 'Score'])
    
    # Load the primer information
    df_primers = pd.read_csv(primer_file)
    
    # Extract serotyping and genotyping from taxonomy
    df_amplicons['Amplicons_Serotyping'] = df_amplicons['Taxonomy'].apply(lambda x: x.split(';')[5])
    df_amplicons['Amplicons_Genotyping'] = df_amplicons['Taxonomy'].apply(lambda x: x.split(';')[6])
    df_whole_genome['Whole_Genome_Serotyping'] = df_whole_genome['Taxonomy'].apply(lambda x: x.split(';')[5])
    df_whole_genome['Whole_Genome_Genotyping'] = df_whole_genome['Taxonomy'].apply(lambda x: x.split(';')[6])
    
    # Initialize list for final data
    final_data = []

    for index, row in df_primers.iterrows():
        primer_pair = row['Primer_Pair']
        amplicon_id = row['Amplicon_ID']
        accession_number = row['Accession_Number']
        sequence = row['Sequence']
        amplicon_length = row['Amplicon_Length']
        
        # Find the amplicons genotyping result
        amplicons_row = df_amplicons[df_amplicons['Amplicons'] == amplicon_id]
        if amplicons_row.empty:
            print(f"Amplicon {amplicon_id} not found in amplicons genotyping DataFrame.")
            continue
        amplicons_serotyping = amplicons_row['Amplicons_Serotyping'].values[0]
        amplicons_genotyping = amplicons_row['Amplicons_Genotyping'].values[0]
        
        # Find the whole genome genotyping result
        whole_genome_row = df_whole_genome[df_whole_genome['Accession_Number'] == accession_number]
        if whole_genome_row.empty:
            print(f"Accession number {accession_number} not found in whole genome genotyping DataFrame.")
            continue
        whole_genome_serotyping = whole_genome_row['Whole_Genome_Serotyping'].values[0]
        whole_genome_genotyping = whole_genome_row['Whole_Genome_Genotyping'].values[0]

        # Append data to the final list
        final_data.append({
            'Primer_Pair': primer_pair,
            'Amplicon_ID': amplicon_id,
            'Accession_Number': accession_number,
            'Sequence': sequence,
            'Amplicon_Length': amplicon_length,
            'Amplicons_Serotyping': amplicons_serotyping,
            'Amplicons_Genotyping': amplicons_genotyping,
            'Whole_Genome_Serotyping': whole_genome_serotyping,
            'Whole_Genome_Genotyping': whole_genome_genotyping
        })

    # Convert the final data list to a DataFrame
    final_df = pd.DataFrame(final_data)

    # Debugging: Print the final DataFrame before saving
    print("Final DataFrame:\n", final_df)
    
    # Save the DataFrame to a CSV file
    final_df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python generate_csv_confusion_matrix.py amplicons_genotyping_file whole_genome_genotyping_file primer_file output_csv")
        sys.exit(1)
    
    amplicons_genotyping_file = sys.argv[1]
    whole_genome_genotyping_file = sys.argv[2]
    primer_file = sys.argv[3]
    output_csv = sys.argv[4]
    generate_final_csv(amplicons_genotyping_file, whole_genome_genotyping_file, primer_file, output_csv)
