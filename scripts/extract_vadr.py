import pandas as pd
import sys

def extract_passed_sequences(sqa_file, sqc_file, original_fasta, output_fasta):
    # Load combined .sqa and .sqc files
    sqa_df = pd.read_csv(sqa_file, delim_whitespace=True, header=None, comment='#', low_memory=False)
    sqc_df = pd.read_csv(sqc_file, delim_whitespace=True, header=None, comment='#', low_memory=False)

    # Assign column names based on VADR documentation
    sqa_df.columns = ["seq_idx_sqa", "seq_name", "seq_len", "p/f", "ant", "model1", "grp1", "sub_grp1", "nfa", "nfn", "nf5", "nf3", "nfalt", "alerts_sqa"]
    sqc_df.columns = ["seq_idx_sqc", "seq_name", "seq_len", "p/f", "ant", "model1", "grp1", "sub_grp1", "score", "sc/nt", "seq_cov", "mdl_cov", "mdl_bias", "mdl_hits", "mdl_str", "model2", "grp2", "sub_grp2", "score_diff", "diff_nt", "seq_alerts_sqc"]

    # Merge the dataframes on sequence name
    combined_df = pd.merge(sqa_df, sqc_df, on="seq_name")

    # Filter out failed sequences
    passed_sequences = combined_df[combined_df["p/f_x"] == "PASS"]["seq_name"].tolist()

    # Read the original FASTA file and write only the passed sequences to the new FASTA file
    with open(original_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        write_flag = False
        for line in infile:
            if line.startswith(">"):
                seq_name = line[1:].strip()
                write_flag = seq_name in passed_sequences
            if write_flag:
                outfile.write(line)

    print(f"Passed sequences extracted to {output_fasta}")

if __name__ == "__main__":
    sqa_file = sys.argv[1]
    sqc_file = sys.argv[2]
    original_fasta = sys.argv[3]
    output_fasta = sys.argv[4]
    extract_passed_sequences(sqa_file, sqc_file, original_fasta, output_fasta)
