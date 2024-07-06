import os
import sys

def split_fasta(file_path, output_dir, max_sequences=20000):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(file_path, 'r') as file:
        seq_count = 0
        file_count = 1
        output_file = os.path.join(output_dir, f'part_{file_count}.fasta')
        out_f = open(output_file, 'w')
        
        for line in file:
            if line.startswith('>'):
                if seq_count >= max_sequences:
                    out_f.close()
                    file_count += 1
                    seq_count = 0
                    output_file = os.path.join(output_dir, f'part_{file_count}.fasta')
                    out_f = open(output_file, 'w')
                seq_count += 1
            out_f.write(line)
        
        out_f.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python split_fasta.py <input_fasta> <output_dir> <max_sequences>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_dir = sys.argv[2]
    max_sequences = int(sys.argv[3])

    split_fasta(input_fasta, output_dir, max_sequences)
