import sys
import logging

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def link_nested_amplicons(outer_results_file, inner_results_file, inner_amplicons_file, primers_file, output_file, min_length=100, max_length=500):
    setup_logging()
    logging.info(f"Linking nested amplicons between {min_length} and {max_length} bp")

    # Read the lengths and accession numbers from the outer results file
    outer_amplicon_data = {}
    with open(outer_results_file, 'r') as res_file:
        header = next(res_file).strip().split('\t')
        if len(header) < 4:
            logging.error("Outer results file header does not have enough columns. Expected 4 columns.")
            sys.exit(1)

        for line in res_file:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                amp_id = parts[0]
                accession_number = parts[1]
                outer_amplicon_data[amp_id] = accession_number
                logging.debug(f"Outer amplicon {amp_id} linked to accession number {accession_number}")
            else:
                logging.warning(f"Skipping malformed line in outer results: {line.strip()}")

    logging.info(f"Total outer amplicon entries read: {len(outer_amplicon_data)}")

    # Read the inner results file and link with the outer results
    inner_amplicon_data = {}
    with open(inner_results_file, 'r') as res_file:
        header = next(res_file).strip().split('\t')
        if len(header) < 4:
            logging.error("Inner results file header does not have enough columns. Expected 4 columns.")
            sys.exit(1)

        for line in res_file:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                amp_id = parts[0]
                outer_amp_id = parts[1]
                try:
                    length = int(parts[3])
                    if outer_amp_id in outer_amplicon_data:
                        accession_number = outer_amplicon_data[outer_amp_id]
                        inner_amplicon_data[amp_id] = (accession_number, length)
                        logging.debug(f"Inner amplicon {amp_id} linked to accession number {accession_number} with length {length}")
                except ValueError:
                    logging.warning(f"Skipping line due to ValueError: {line.strip()}")
            else:
                logging.warning(f"Skipping malformed line in inner results: {line.strip()}")

    logging.info(f"Total inner amplicon entries read: {len(inner_amplicon_data)}")

    # Read the primer pairs information
    primer_data = {}
    with open(primers_file, 'r') as primers_file:
        for line in primers_file:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                primer_pair = parts[2]
                primer_data[primer_pair] = parts[:2]
                logging.debug(f"Read primer pair {primer_pair} with parts: {parts[:2]}")

    logging.info(f"Total primer pair entries read: {len(primer_data)}")

    logging.info("Starting to link amplicons")
    # Filter the amplicons based on the lengths and include accession numbers and primer pairs in the output
    with open(inner_amplicons_file, 'r') as amp_file, open(output_file, 'w') as out_file:
        for line in amp_file:
            if line.startswith('>'):
                amp_id = line[1:].strip()
                if amp_id in inner_amplicon_data:
                    accession_number, length = inner_amplicon_data[amp_id]
                    if min_length <= length <= max_length:
                        for primer_pair, primers in primer_data.items():
                            if amp_id.startswith(primer_pair):
                                out_file.write(f"{primer_pair},{amp_id},{accession_number}\n")
                                logging.debug(f"Writing: {primer_pair},{amp_id},{accession_number}")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python link_nested_amplicons.py outer_results_file inner_results_file inner_amplicons_file primers_file output_file")
        sys.exit(1)
    setup_logging()
    outer_results_file = sys.argv[1]
    inner_results_file = sys.argv[2]
    inner_amplicons_file = sys.argv[3]
    primers_file = sys.argv[4]
    output_file = sys.argv[5]
    link_nested_amplicons(outer_results_file, inner_results_file, inner_amplicons_file, primers_file, output_file)
