import sys
import re
import os
import argparse
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description="Assign taxonomy based on BLAST results.")
    parser.add_argument("blast_files", nargs='+', help="BLAST output files.")
    parser.add_argument("tax_file", help="Taxonomy file.")
    parser.add_argument("output_file", help="Output file for taxonomy assignment.")
    return parser.parse_args()

def parse_taxonomy_file(tax_file):
    taxonomy_dict = {}
    with open(tax_file, 'r') as file:
        for line in file:
            elements = line.strip().split('\t')
            if len(elements) != 2:
                print(f"Skipping invalid line in taxonomy file: {line.strip()}")
                continue
            tax_id, tax_info = elements
            taxonomy_dict[tax_id] = re.sub(r'D_[0-9]*__', '', tax_info).replace(' ', '_').split(';')
    return taxonomy_dict

def parse_blast_file(blast_file):
    blast_hits = defaultdict(list)
    with open(blast_file, 'r') as file:
        for line in file:
            elements = line.strip().split('\t')
            if len(elements) < 12:
                continue
            qseqid, qlen, sseqid, pident, length, qstart, qend, sstart, send, evalue, bitscore, staxids = elements
            blast_hits[qseqid].append((float(pident), sseqid, staxids))
    return blast_hits

def assign_taxonomy(blast_hits, taxonomy_dict):
    assigned_taxonomy = {}
    for query, hits in blast_hits.items():
        if not hits:
            assigned_taxonomy[query] = ['unassigned'] * 7
            continue
        best_hit = max(hits, key=lambda x: x[0])
        _, best_sseqid, best_taxid = best_hit
        if best_taxid in taxonomy_dict:
            assigned_taxonomy[query] = taxonomy_dict[best_taxid]
        else:
            print(f"No match found in taxonomy dictionary for taxid: {best_taxid}")
            assigned_taxonomy[query] = ['unassigned'] * 7
    return assigned_taxonomy

def write_output(assigned_taxonomy, output_file):
    with open(output_file, 'w') as file:
        for query, taxonomy in assigned_taxonomy.items():
            file.write(f"{query}\t{';'.join(taxonomy)}\n")

def main():
    args = parse_arguments()
    print("###Program Start###\n")
    
    print("Parsing taxonomy file...")
    taxonomy_dict = parse_taxonomy_file(args.tax_file)
    print("Parsed taxonomy entries:", len(taxonomy_dict))
    
    print("Parsing BLAST results and assigning taxonomy...")
    all_blast_hits = defaultdict(list)
    for blast_file in args.blast_files:
        blast_hits = parse_blast_file(blast_file)
        for query, hits in blast_hits.items():
            all_blast_hits[query].extend(hits)
    
    print("Total queries in BLAST results:", len(all_blast_hits))
    assigned_taxonomy = assign_taxonomy(all_blast_hits, taxonomy_dict)
    
    print("Writing output file...")
    write_output(assigned_taxonomy, args.output_file)
    
    print("###Program End###\n")

if __name__ == "__main__":
    main()
