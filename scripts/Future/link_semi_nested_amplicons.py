import sys
import logging

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def read_results_file(results_file):
    data = {}
    with open(results_file, 'r') as file:
        next(file)  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                amp_id = parts[0]
                accession_number = parts[1]
                try:
                    length = int(parts[3])
                    data[amp_id] = (accession_number, length)
                except ValueError:
                    logging.warning(f"Skipping line due to ValueError: {line.strip()}")
    return data

def link_amplicons(outer_results_file, inner_results_file, amplicons_file, primers_file, output_file, min_length=100, max_length=500):
    setup_logging()
    
    # Read the outer and inner results
    logging.info("Reading outer results file")
    outer_data = read_results_file(outer_results_file)
    
    logging.info("Reading inner results file")
    inner_data = read_results_file(inner_results_file)
    
    # Read the primer pairs information
    logging.info("Reading primer pairs information")
    primer_data = {}
    with open(primers_file, 'r') as primers_file:
        for line in primers_file:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                primer_pair = parts[2]
                primer_data[primer_pair] = parts[:2]

    # Link the amplicons
    logging.info("Linking amplicons")
    with open(amplicons_file, 'r') as amp_file, open(output_file, 'w') as out_file:
        write = False
        for line in amp_file:
            if line.startswith('>'):
                write = False
                amp_id = line[1:].strip()
                if amp_id in inner_data:
                    inner_acc_num, length = inner_data[amp_id]
                    if inner_acc_num in outer_data:
                        accession_number = outer_data[inner_acc_num][0]
                        if min_length <= length <= max_length:
                            write = True
                            out_file.write(f"{amp_id},{inner_acc_num},{accession_number}\n")
                            logging.debug(f"Writing: {amp_id},{inner_acc_num},{accession_number}")
                        else:
                            logging.debug(f"Skipping amplicon {amp_id} with length {length} (out of range)")
                    else:
                        logging.debug(f"Skipping amplicon {amp_id} (not found in outer results)")
                else:
                    logging.debug(f"Skipping amplicon {amp_id} (not found in inner results)")
            elif write:
                # Do not write the sequence part
                write = False

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python link_semi_nested_amplicons.py outer_results_file inner_results_file amplicons_file primers_file output_file")
        sys.exit(1)
    
    outer_results_file = sys.argv[1]
    inner_results_file = sys.argv[2]
    amplicons_file = sys.argv[3]
    primers_file = sys.argv[4]
    output_file = sys.argv[5]
    
    link_amplicons(outer_results_file, inner_results_file, amplicons_file, primers_file, output_file)
import sys
import logging

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def read_results_file(results_file):
    data = {}
    with open(results_file, 'r') as file:
        next(file)  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                amp_id = parts[0]
                accession_number = parts[1]
                try:
                    length = int(parts[3])
                    data[amp_id] = (accession_number, length)
                except ValueError:
                    logging.warning(f"Skipping line due to ValueError: {line.strip()}")
    return data

def link_amplicons(outer_results_file, inner_results_file, amplicons_file, primers_file, output_file, min_length=100, max_length=500):
    setup_logging()
    
    # Read the outer and inner results
    logging.info("Reading outer results file")
    outer_data = read_results_file(outer_results_file)
    
    logging.info("Reading inner results file")
    inner_data = read_results_file(inner_results_file)
    
    # Read the primer pairs information
    logging.info("Reading primer pairs information")
    primer_data = {}
    with open(primers_file, 'r') as primers_file:
        for line in primers_file:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                primer_pair = parts[2]
                primer_data[primer_pair] = parts[:2]

    # Link the amplicons
    logging.info("Linking amplicons")
    with open(amplicons_file, 'r') as amp_file, open(output_file, 'w') as out_file:
        write = False
        for line in amp_file:
            if line.startswith('>'):
                write = False
                amp_id = line[1:].strip()
                if amp_id in inner_data:
                    inner_acc_num, length = inner_data[amp_id]
                    if inner_acc_num in outer_data:
                        accession_number = outer_data[inner_acc_num][0]
                        if min_length <= length <= max_length:
                            write = True
                            out_file.write(f"{amp_id},{inner_acc_num},{accession_number}\n")
                            logging.debug(f"Writing: {amp_id},{inner_acc_num},{accession_number}")
                        else:
                            logging.debug(f"Skipping amplicon {amp_id} with length {length} (out of range)")
                    else:
                        logging.debug(f"Skipping amplicon {amp_id} (not found in outer results)")
                else:
                    logging.debug(f"Skipping amplicon {amp_id} (not found in inner results)")
            elif write:
                # Do not write the sequence part
                write = False

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python link_semi_nested_amplicons.py outer_results_file inner_results_file amplicons_file primers_file output_file")
        sys.exit(1)
    
    outer_results_file = sys.argv[1]
    inner_results_file = sys.argv[2]
    amplicons_file = sys.argv[3]
    primers_file = sys.argv[4]
    output_file = sys.argv[5]
    
    link_amplicons(outer_results_file, inner_results_file, amplicons_file, primers_file, output_file)
