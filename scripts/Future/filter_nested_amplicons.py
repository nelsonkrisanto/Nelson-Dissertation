import sys
import logging

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def filter_amplicons(results_file, amplicons_file, output_file, min_length=100, max_length=500):
    setup_logging()
    logging.info(f"Filtering amplicons between {min_length} and {max_length} bp")

    # Read the lengths from the results file
    lengths = {}
    with open(results_file, 'r') as res_file:
        header = next(res_file).strip().split('\t')
        if len(header) < 4:
            logging.error("Results file header does not have enough columns. Expected 4 columns.")
            sys.exit(1)

        for line in res_file:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                logging.warning(f"Skipping malformed line: {line.strip()}")
                continue

            amp_id, length_str = parts[0], parts[3]
            try:
                length = int(length_str)
                lengths[amp_id] = length
                logging.debug(f"Read length {length} for amplicon {amp_id}")
            except ValueError:
                logging.warning(f"Skipping line due to ValueError: {line.strip()}")

    # Filter the amplicons based on the lengths
    logging.info("Starting to filter amplicons")
    with open(amplicons_file, 'r') as amp_file, open(output_file, 'w') as out_file:
        write = False
        for line in amp_file:
            if line.startswith('>'):
                write = False
                amp_id = line[1:].strip()
                if amp_id in lengths:
                    length = lengths[amp_id]
                    if min_length <= length <= max_length:
                        write = True
                        out_file.write(line)
                        logging.debug(f"Writing amplicon {amp_id} with length {length}")
                    else:
                        logging.debug(f"Skipping amplicon {amp_id} with length {length} (out of range)")
                else:
                    logging.debug(f"Skipping amplicon {amp_id} (not found in results)")
            elif write:
                out_file.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_amplicons.py results_file amplicons_file output_file")
        sys.exit(1)
    results_file = sys.argv[1]
    amplicons_file = sys.argv[2]
    output_file = sys.argv[3]
    filter_amplicons(results_file, amplicons_file, output_file)
