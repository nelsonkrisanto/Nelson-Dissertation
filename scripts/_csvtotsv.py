import csv

def convert_csv_to_tsv(input_csv_path, output_tsv_path):
    # Open the input CSV file
    with open(input_csv_path, newline='') as csv_file:
        # Create a CSV reader object specifying the delimiter as a comma
        csv_reader = csv.reader(csv_file, delimiter=',')
        
        # Open the output TSV file
        with open(output_tsv_path, 'w', newline='') as tsv_file:
            # Create a CSV writer object specifying the delimiter as a tab
            tsv_writer = csv.writer(tsv_file, delimiter='\t')
            
            # Iterate through the CSV file and write each row to the TSV file
            for row in csv_reader:
                tsv_writer.writerow(row)

# Specify the full or relative path to your input CSV file
input_csv_path = 'C:/Users/Nelso/OneDrive/Documents/Thesis/data/primer_metadata.csv'

# Specify the full or relative path to your desired output TSV file
output_tsv_path = 'C:/Users/Nelso/OneDrive/Documents/Thesis/data/primer_metadata.tsv'

# Call the function with the specified paths
convert_csv_to_tsv(input_csv_path, output_tsv_path)

print("Done Converting...")