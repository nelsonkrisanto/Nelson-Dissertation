import os

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

# Define the paths
fasta_dir="/home/people/23203786/scratch/Nelson-Dissertation/raw_data/NCBI_dengue_data/"
fasta_file = os.path.join(fasta_dir, "1.1.dengue_virus_sequences.fasta")
output_dir = "/home/people/23203786/scratch/Nelson-Dissertation/results/fasta"

# Run the split function
split_fasta(fasta_file, output_dir)
