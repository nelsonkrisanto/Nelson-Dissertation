import sys

def filter_amplicons(results_file, amplicons_file, output_file, min_length=300, max_length=1200):
    # Read the lengths from the results file
    lengths = {}
    with open(results_file, 'r') as res_file:
        next(res_file)  # Skip the header line
        for line in res_file:
            parts = line.strip().split('\t')
            if len(parts) >= 4:  # Ensure there are enough parts in the line
                amp_id = parts[0]
                try:
                    length = int(parts[3])
                    lengths[amp_id] = length
                except ValueError:
                    print(f"Skipping line due to ValueError: {line.strip()}")
            else:
                print(f"Skipping malformed line: {line.strip()}")

    # Filter the amplicons based on the lengths
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
