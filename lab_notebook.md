# Title
## High-Level Guide to Using the Scripts

### Configuration: `00.config.sh`
Please update the script with your specific settings where you see the comment "## adjust". Ensure these settings align with your project requirements.

### Fetching Dengue Virus Sequence Data: `1.1.fetch_NCBI.sh`
Forked from 'Benchmarking_PCR_ONT' by Amy Fitzpatrick. The script has been modified as follows:
- Modify the sections marked as "## adjust" to suit your research needs.
- When you reach the portion of the script that says, "# Perform the search to get dengue virus sequences. Adjust the search term as needed," refer to the instructions provided below.

To tailor the search query for your need, use the NCBI Nucleotide Advanced Search Builder. The image below (Figure 1.1) illustrates the Advanced Search Builder interface, which you can access at the following URL: [NCBI Nucleotide Advanced Search Builder](https://www.ncbi.nlm.nih.gov/nuccore/advanced).

![Figure 1.1 - NCBI Nucleotide Advanced Search Builder](images/01_01.png)
<p align="center">
  Figure 1.1 - NCBI Nucleotide Advanced Search Builder
</p>

In the Advanced Search Builder, adjust the search parameters as necessary. Your configured parameters will reflect in the search details field, as highlighted in Figure 1.1.
- **Output File**:  
  - `1.1.dengue_sequence_info.tsv`
  - `1.1.dengue_virus_sequences.fasta`

### Fetching Reference Genome: `2.1.fetch_ref_genome.sh`
Retrieve high-quality reference genomes of the Dengue virus from the NCBI nucleotide database.
- Access NCBI nucleotide database targeting RefSeq entries.
- Ensure compliance with ethical data fetching standards.
- Store downloaded genomes in FASTA format in a designated directory
- **Output File**: 
  - `GCF_000862125.1_ViralProj15306_genomic.fna`
  - `GCF_000865065.1_ViralProj15599_genomic.fna`
  - `GCF_000866625.1_ViralProj15598_genomic.fna`
  - `GCF_000871845.1_ViralProj20183_genomic.fna`

### Indexing Reference Genome: `2.2.index_ref_genome.sh`
Index the retrieved Dengue virus reference genomes to facilitate efficient sequence alignments.
- **Input File**: 
  - `GCF_000862125.1_ViralProj15306_genomic.fna`
  - `GCF_000865065.1_ViralProj15599_genomic.fna`
  - `GCF_000866625.1_ViralProj15598_genomic.fna`
  - `GCF_000871845.1_ViralProj20183_genomic.fna`
- Sequentially handle each genome file.
- Use the `bwa index` command to create indices.
- Store indexing results in the same directory as the genome files.

### Fetching NCBI Data: `3.1.scrappaper.py`
Forked from 'ScrapPaper' by M. R. Rafsanjani. The script has been modified as follows:
- Removed the Google Scholar function.
- Added a function to scrape from both NCBI Pubmed Summary and Abstract formats.
- Added a function to merge the results.
- Added a function to extract the year from references.
- Modified the pandas write section to append new data instead of replacing it.

### Cleanning PCR Primers Sequence from the NCBI Data Result
This requires manual work. Go through the results from the '3.1.scrappaper.py' section and look for the Dengue PCR primers. Follow these guidelines:
- **Input File**: `all_cleanned_v1.xlsx`
- Only Dengue PCR Primers will be taken.
- Remove publications with no primers or where the result is N/A.
- Rename primers with no name to the first letter of every word in the publication name followed by a number.
- Remove artificial tags, probes, and RT elements.
- Replace spaces with underscores for primer names.
- Remove primer with no known genotype
- **Output File**: `cleaned_primers.csv`

### Convert Cleaned CSV PCR Primers Result to FASTA Format: `3.2.csvtofasta.py`
- **Input File**: `cleaned_primers.csv`
- **Output File**: `dengue_primers.fasta`

### Check Primer Position: `3.3.check_primer_position.sh`
Forked from 'Benchmarking_PCR_ONT' by Amy Fitzpatrick. The script has been modified as follows:
- Added several error checks to ensure the accuracy and integrity of the primer positions.

- **Input File**: `dengue_primers.fasta`
- **Output File**: `mapping_positions.tsv`
- **Functionality**:
  - Load primer sequences from a FASTA file.
  - Align each primer to reference genomes using BWA.
  - Extract alignment positions and save to a TSV file.
  - Added error checks to ensure the accuracy of the alignments.

### Fetch Temperature: `3.4.fetch_temp.sh`
This script calculates the minimum and maximum melting temperatures (Tm) for each primer sequence using Wallace's rule and the nearest neighbor method from an Excel file containing primer sequences. The results are then saved in a TSV file.
- **Input File**: `cleanned_primers.xlsx`
- **Output File**: `primer_metadata.tsv`
- **Functionality**:
  - Load primer sequences from an Excel file.
  - Handle ambiguous bases in primer sequences.
  - Generate all possible sequences for ambiguous primers.
  - Calculate average minimum and maximum melting temperatures (Tm) for each sequence.
  - Save the results to a TSV file.

### Generation of Primers Pairs

### Conventional PCR Primer Combination: `3.5.1.conventional.sh`
- **Input**: `mapping_positions.tsv`, `primer_metadata.tsv`
- **Output**: `conventional_primer_combinations.tsv`
- The script generates primer combinations for conventional PCR by pairing forward and reverse primers based on specified criteria.

### Nested PCR Primer Combination: `3.5.2.nested.sh`
- **Input**: `mapping_positions.tsv`, `primer_metadata.tsv`
- **Output**: `all_nested_primer_combinations.tsv`, `top_200_nested_primer_combinations.tsv`, `matched_nested_primer_combinations.tsv`
- This script generates nested PCR primer combinations in two rounds:
  1. First round: Identifies valid primer pairs for the outer PCR reaction.
  2. Second round: Identifies valid primer pairs for the inner PCR reaction that fall within the outer amplicon.
- **Parameters and Requirements:**
  - Minimum primer distance: 50 bases.
  - Melting temperature difference (Tm): ≤ 5°C.
  - GC content difference: ≤ 10%.
  - Homopolymer runs: ≤ 4 bases.
  - Valid primer pairs must be within the specified regions of interest.
- **Output Files:**
  - `all_nested_primer_combinations.tsv`: All valid nested primer combinations.
  - `top_200_nested_primer_combinations.tsv`: Top 200 nested primer combinations based on specified criteria.
  - `matched_nested_primer_combinations.tsv`: Matched combinations of first and second round primers.

### Semi-Nested PCR Primer Combination: `3.5.3.semi_nested.sh`
- **Input**: `mapping_positions.tsv`, `primer_metadata.tsv`
- **Output**: `all_semi_nested_primer_combinations_strict.tsv`, `top_200_semi_nested_primer_combinations_strict.tsv`, `matched_semi_nested_primer_combinations_strict.tsv`
- This script generates semi-nested PCR primer combinations in two rounds:
  1. First round: Identifies valid primer pairs for the outer PCR reaction.
  2. Second round: Identifies valid new forward or reverse primers for the semi-nested PCR that fall within the outer amplicon.
- **Parameters and Requirements:**
  - Minimum primer distance: 50 bases.
  - Melting temperature difference (Tm): ≤ 3°C (stricter).
  - GC content difference: ≤ 10%.
  - Homopolymer runs: ≤ 4 bases.
  - Valid primer pairs must be within the specified regions of interest.
- **Output Files:**
  - `all_semi_nested_primer_combinations_strict.tsv`: All valid semi-nested primer combinations.
  - `top_200_semi_nested_primer_combinations_strict.tsv`: Top 200 semi-nested primer combinations based on specified criteria.
  - `matched_semi_nested_primer_combinations_strict.tsv`: Matched combinations of first and second round primers.

### Phylogenetic Tree Creation: `phylogenetic_tree.ipynb`
- Download nucleotide sequences from NCBI for a list of the following accession numbers:
  - DQ235145
  - AY323490
  - AF331718
  - AF253419
  - Y07863
  - D12937
  - X86784
  - DQ235152
  - DQ235151
  - DQ235153
  - AY193805
  - L06436
  - AF311056
  - DQ235149
  - U27495
  - X07755
  - L40361
  - DQ235144
  - DQ235150
  - KF815939
  - DQ235146
  - AY632536
  - AF013366
  - AF013375
  - AF013390
  - U88536
  - U87411
  - M93130
  - AF326573
  - KF917536
  - M18370
  - AF013384
  - AF161266
  - DQ525916
  - AY453411
  - D00246
  - M12294
  - AF013413
  - AY632541
  - AF013407
  - AY632545
  - AY632539
  - AF013397
  - AF013377
  - JX236040
  - JF895923
  - AY632535
  - DQ837642
  - EU707555
  - X03700
  - AY632540
  - DQ859056
  - DQ859057
  - DQ859060
  - DQ859066
  - DQ859067
  - DQ859062
  - DQ859065
  - DQ837641
  - AF013405
  - AB114858
  - AF160193
  - AF013370
  - KJ469371
  - AJ242984
  - AF013401
  - AF013402
  - AF013365
  - AF013368
  - AF013371
  - AJ299445
  - AF013369
  - AF013394
  - AF144692
- **Manual Work**:
  - Replace spaces with underscores in result.
  - Remove empty lines.
  - Use [Clustal Omega](https://www.ebi.ac.uk/jdispatcher/msa/muscle/) for alignment.
  - Create the phylogenetic tree using iTOL.
- **Output Files:**: 
  - `orthoflavivirus_sequences.fasta`
  - `tree.......... to be updated`

### In Silico PCR
- The script simulates PCR amplification using specified primers to validate their performance.

### Pre-requisite/Dependencies
The following modules, environment, and packages are needed in this project:
* **Modules:**
  - Anaconda
* **Environment:**
  - Conda Environment
* **Packages:**
  - Entrez Direct (EDirect)
  - Python libraries: pandas, BeautifulSoup, requests, Biopython, tqdm, logging

By following this structured approach, you can effectively generate, validate, and analyze PCR primers for Dengue virus research using both nested and semi-nested PCR methods.
