import pandas as pd
import sys

# List of other orthoflavivirus accession numbers
other_orthoflavivirus_accessions = [
    "DQ235145", "AY323490", "AF331718", "AF253419", "Y07863",
    "D12937", "X86784", "DQ235152", "DQ235151", "DQ235153",
    "AY193805", "L06436", "AF311056", "DQ235149", "U27495",
    "X07755", "L40361", "DQ235144", "DQ235150", "KF815939",
    "DQ235146", "AY632536", "AF013366", "AF013375", "AF013390",
    "U88536", "U87411", "M93130", "AF326573", "KF917536",
    "M18370", "AF013384", "AF161266", "DQ525916", "AY453411",
    "D00246", "M12294", "AF013413", "AY632541", "AF013407",
    "AY632545", "AY632539", "AF013397", "AF013377", "JX236040",
    "JF895923", "AY632535", "DQ837642", "EU707555", "X03700",
    "AY632540", "DQ859056", "DQ859057", "DQ859060", "DQ859066",
    "DQ859067", "DQ859062", "DQ859065", "DQ837641", "AF013405",
    "AB114858", "AF160193", "AF013370", "KJ469371", "AJ242984",
    "AF013401", "AF013402", "AF013365", "AF013368", "AF013371",
    "AJ299445", "AF013369", "AF013394", "AF144692"
]

def generate_final_csv(amplicons_genotyping_file, whole_genome_genotyping_file, primer_file, output_csv):
    # Load the taxonomy assignment results for amplicons and whole genome
    df_amplicons = pd.read_csv(amplicons_genotyping_file, sep='\t', header=None, names=[
        'Amplicons', 'Taxonomy', 'Score'])
    df_whole_genome = pd.read_csv(whole_genome_genotyping_file, sep='\t', header=None, names=[
        'Accession_Number', 'Taxonomy', 'Score'])
    
    # Load the primer information
    df_primers = pd.read_csv(primer_file, header=None, names=[
        'Primer_Pairs', 'Amplicons', 'Accession_Number'])
    
    # Extract genotyping from taxonomy
    df_amplicons['Amplicons_Genotyping'] = df_amplicons['Taxonomy'].apply(lambda x: x.split(';')[5])
    df_whole_genome['Whole_Genome_Genotyping'] = df_whole_genome['Taxonomy'].apply(lambda x: x.split(';')[5])
    
    # Initialize list for final data
    final_data = []

    for index, row in df_primers.iterrows():
        primer_pair = row['Primer_Pairs']
        amplicon = row['Amplicons']
        accession_number = row['Accession_Number']
        
        # Find the amplicons genotyping result
        amplicons_row = df_amplicons[df_amplicons['Amplicons'] == amplicon]
        if amplicons_row.empty:
            print(f"Amplicon {amplicon} not found in amplicons genotyping DataFrame.")
            continue
        amplicons_genotyping = amplicons_row['Amplicons_Genotyping'].values[0]
        
        # Find the whole genome genotyping result
        whole_genome_row = df_whole_genome[df_whole_genome['Accession_Number'] == accession_number]
        if whole_genome_row.empty:
            print(f"Accession number {accession_number} not found in whole genome genotyping DataFrame.")
            continue
        whole_genome_genotyping = whole_genome_row['Whole_Genome_Genotyping'].values[0]

        # Append data to the final list
        final_data.append({
            'Primer Pairs': primer_pair,
            'Amplicons': amplicon,
            'Accession Number': accession_number,
            'Amplicons Genotyping': amplicons_genotyping,
            'Whole Genome Genotyping': whole_genome_genotyping
        })

    # Add dummy entries for other orthoflavivirus sequences
    for accession in other_orthoflavivirus_accessions:
        final_data.append({
            'Primer Pairs': 'dummy_primer',
            'Amplicons': accession,
            'Accession Number': accession,
            'Amplicons Genotyping': 'other_virus',
            'Whole Genome Genotyping': 'other_virus'
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
