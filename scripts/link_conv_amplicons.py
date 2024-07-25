import os
import sys

def link_conv_amplicons(results_file, amplicons_file, primers_file, output_file, min_length=300, max_length=1200):
    # Read the lengths and accession numbers from the results file
    amplicon_data = {}
    with open(results_file, 'r') as res_file:
        next(res_file)  # Skip the header line
        for line in res_file:
            parts = line.strip().split('\t')
            if len(parts) >= 4:  # Ensure there are enough parts in the line
                amp_id = parts[0]
                accession_number = parts[1]
                try:
                    length = int(parts[3])
                    amplicon_data[amp_id] = (accession_number, length)
                except ValueError:
                    # Skip lines that don't have a valid integer length
                    continue
            else:
                # Skip lines that do not have the correct number of parts
                continue

    # Read the primer pairs information
    primer_data = {}
    with open(primers_file, 'r') as primers_file:
        for line in primers_file:
            parts = line.strip().split('\t')
            if len(parts) >= 3:  # Ensure there are enough parts in the line
                primer_pair = parts[2]
                primer_data[primer_pair] = parts[:2]

    # Read the sequences from the amplicons file
    sequence_data = {}
    with open(amplicons_file, 'r') as amp_file:
        sequence_id = None
        sequence = ""
        for line in amp_file:
            if line.startswith('>'):
                if sequence_id:
                    sequence_data[sequence_id] = sequence
                sequence_id = line[1:].strip()
                sequence = ""
            else:
                sequence += line.strip()
        if sequence_id:
            sequence_data[sequence_id] = sequence

    # Filter the amplicons based on the lengths and include accession numbers and primer pairs in the output
    with open(output_file, 'w') as out_file:
        out_file.write("Primer_Pair,Amplicon_ID,Accession_Number,Sequence,Amplicon_Length\n")
        for amp_id, (accession_number, length) in amplicon_data.items():
            if min_length <= length <= max_length:
                for primer_pair, primers in primer_data.items():
                    if amp_id.startswith(primer_pair):
                        sequence = sequence_data.get(amp_id, "")
                        out_file.write(f"{primer_pair},{amp_id},{accession_number},{sequence},{length}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python link_conv_amplicons.py results_file amplicons_file primers_file output_file")
        sys.exit(1)
    results_file = sys.argv[1]
    amplicons_file = sys.argv[2]
    primers_file = sys.argv[3]
    output_file = sys.argv[4]
    link_conv_amplicons(results_file, amplicons_file, primers_file, output_file)

