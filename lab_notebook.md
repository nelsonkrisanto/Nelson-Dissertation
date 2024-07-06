# Comparison of Bioinformatic Approaches for Minimizing Sequence Redundancy in Amplicon Sequencing Data from Oxford Nanopore Technology Platforms Applied to Dengue Virus Research

## High-Level Guide to Using the Scripts

### Configuration: `00.config.sh`
Ensure you update the script with your specific settings where you see the comment "## adjust". These settings must align with your project requirements.

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

- **Script**: `1.1.fetch_NCBI.sh`
- **Packages and Libraries**: Entrez Direct (version 13.2)
- **Input File**: NCBI nucleotide database query results
- **Output File**: `1.1.dengue_virus_sequences.fasta`

### Clustering Dengue Virus Sequences: `1.2.cluster_MMseqs2.sh`
The script `1.2.cluster_MMseqs2.sh` was used to cluster Dengue virus sequences using MMseqs2 (version 13-45111). The script created a database from the input FASTA file, performed clustering with a 97% identity threshold, and generated an output FASTA file containing representative sequences of each cluster.

- **Script**: `1.2.cluster_MMseqs2.sh`
- **Packages and Libraries**: MMseqs2 (version 13-45111)
- **Input File**: `1.1.dengue_virus_sequences.fasta`
- **Output Files**:
  - MMseqs2 database: `db/input_db`
  - Clustering results: `clusters`
  - Clustered sequences FASTA file: `clustered_sequences.fasta`

### Generating VADR Amplicons: `1.3.vadr_amplicons.sh`
The script `1.3.vadr_amplicons.sh` generated amplicons for VADR analysis. This script set up the environment for Perl scripts required for the analysis and utilized the `v-annotate.pl` script using VADR models specific to Flaviviridae from the VADR library (https://github.com/ncbi/vadr).

- **Script**: `1.3.vadr_amplicons.sh`
- **Packages and Libraries**: Anaconda (version 3.5.2), VADR (version 1.5.1), Perl
- **Input File**: `filtered_conv_insilico_PCR_amplicons.fasta`
- **Output Files**:
  - VADR output directory: `vadr_output`
  - VADR annotated files

### Extracting VADR Amplicons: `1.4.extract_vadr_amplicons.sh`
Following the generation of amplicons, the script `1.4.extract_vadr_amplicons.sh` was employed to extract VADR amplicons. This script merged all `.sqa` and `.sqc` files and used a Python script (`extract_vadr.py`) to filter and retain only the sequences that passed quality checks.

- **Script**: `1.4.extract_vadr_amplicons.sh`
- **Packages and Libraries**: Python (version 3.9.15), Perl
- **Tools**: `cat` (for merging files), custom Python script `extract_vadr.py`
- **Input Files**:
  - Merged `.sqa` files: `vadr_combined.sqa`
  - Merged `.sqc` files: `vadr_combined.sqc`
  - Original FASTA file: `filtered_conv_insilico_PCR_amplicons.fasta`
- **Output File**: `passed_sequences.fasta`

### Creating Taxonomy Database: `1.5.taxonomy_db.sh`
The `1.5.taxonomy_db.sh` script is used to create the `protein_db_taxonomy.tsv` file and to combine FASTA files of each genotype from their respective serotypes. The taxonomy data is sourced from NCBI, formatted as follows: Amarillovirales;Flaviviridae;Orthoflavivirus;Orthoflavivirus_denguei;dengue_virus;dengue_virus_type_1;Strain_name. This taxonomy is then mapped to the latest naming convention described in the paper "A new lineage nomenclature to aid genomic surveillance of dengue virus." Additionally, the nomenclature used in Genome Detective is mapped to the new lineage nomenclature, as it will be utilized in the "Assigning Taxonomy Using Blast" step.

#### Nomenclature of Dengue Virus Genotype

| NCBI Taxonomy Entry | New Lineage Nomenclature | Genome Detective Genotype |
|---------------------|--------------------------|---------------------------|
| Dengue virus 1 Brazil/97-11/1997 | DENV-1 genotype V | DENV1_V |
| Dengue virus 1 Jamaica/CV1636/1977 | DENV-1 genotype II | DENV1_II |
| Dengue virus 1 Nauru/West Pac/1974 | DENV-1 genotype I | DENV1_I |
| Dengue virus 1 Singapore/S275/1990 | DENV-1 genotype IV | DENV1_IV |
| Dengue virus 1 Thailand/AHF 82-80/1980 | DENV-1 genotype III | DENV1_III |
| Dengue virus 2 16681-PDK53 | DENV-2 genotype Asian I | DENV2_I |
| Dengue virus 2 China/D2-04 | DENV-2 genotype Asian II | DENV2_II |
| Dengue virus 2 Jamaica/1409/1983 | DENV-2 genotype American | DENV2_AM |
| Dengue virus 2 Malaysia M2 | DENV-2 genotype Asian I | DENV2_I |
| Dengue virus 2 Malaysia M3 | DENV-2 genotype Asian II | DENV2_II |
| Dengue virus 2 Peru/IQT2913/1996 | DENV-2 genotype Cosmopolitan | DENV2_CO |
| Dengue virus 2 Puerto Rico/PR159-S1/1969 | DENV-2 genotype American | DENV2_AM |
| Dengue virus 2 Thailand/0168/1979 | DENV-2 genotype Asian I | DENV2_I |
| Dengue virus 2 Thailand/16681/84 | DENV-2 genotype Asian II | DENV2_II |
| Dengue virus 2 Thailand/NGS-C/1944 | DENV-2 genotype Cosmopolitan | DENV2_CO |
| Dengue virus 2 Thailand/PUO-218/1980 | DENV-2 genotype Asian I | DENV2_I |
| Dengue virus 2 Thailand/TH-36/1958 | DENV-2 genotype Asian II | DENV2_II |
| Dengue virus 2 Tonga/EKB194/1974 | DENV-2 genotype Cosmopolitan | DENV2_CO |
| Dengue virus 3 China/80-2/1980 | DENV-3 genotype III | DENV3_III |
| Dengue virus 3 Martinique/1243/1999 | DENV-3 genotype IV | DENV3_IV |
| Dengue virus 3 Philippines/H87/1956 | DENV-3 genotype I | DENV3_I |
| Dengue virus 3 Singapore/8120/1995 | DENV-3 genotype II | DENV3_II |
| Dengue virus 3 Sri Lanka/1266/2000 | DENV-3 genotype V | DENV3_V |
| Dengue virus 4 Dominica/814669/1981 | DENV-4 genotype II | DENV4_II |
| Dengue virus 4 Philippines/H241/1956 | DENV-4 genotype I | DENV4_I |
| Dengue virus 4 Singapore/8976/1995 | DENV-4 genotype II | DENV4_II |
| Dengue virus 4 Thailand/0348/1991 | DENV-4 genotype I | DENV4_I |
| Dengue virus 4 Thailand/0476/1997 | DENV-4 genotype I | DENV4_I |

- **Script**: `1.5.taxonomy_db.sh`
- **Packages and Libraries**: Python (version 3.7.4)
- **Input Files**: FASTA files
- **Output Files**:
  - Combined FASTA files for each serotype
  - `protein_db_taxonomy.tsv`

### Converting and Creating BLAST Database: `1.6.convert_create_blastdb.sh`
The script `1.6.convert_create_blastdb.sh` converted the combined FASTA files into a format suitable for creating a BLAST database. This process utilized the `fasta2RDP.py` script from the FROGS library (https://github.com/geraldinepascal/FROGS-wrappers) and created a BLAST database.

- **Script**: `1.6.convert_create_blastdb.sh`
- **Packages and Libraries**: Python (version 3.7.4), FROGS library
- **Tools**: `seqkit` (for header manipulation), BLAST (version 2.12.0)
- **Input Files**:
  - Combined FASTA file: `combined_denv4.fasta`
  - Taxonomy file: `protein_db_taxonomy.tsv`
- **Output Files**:
  - RDP formatted FASTA and taxonomy files
  - BLAST database

### Assigning Taxonomy Using BLAST: `1.7.blast_assign_taxonomy.sh`
The script `1.7.blast_assign_taxonomy.sh` assigned taxonomy to the sequences using BLAST. This script ran BLAST searches and used the `taxonomy_assignment_blast.py` script from Joseph7e's repository (https://github.com/Joseph7e/Assign-Taxonomy-with-BLAST) to assign taxonomy based on BLAST results.

- **Script**: `1.7.blast_assign_taxonomy.sh`
- **Packages and Libraries**: BLAST (version 2.12.0), Biopython
- **Tools**: `seqkit` (for removing duplicate sequences)
- **Input Files**:
  - Query FASTA file: `passed_sequences.fasta`
  - BLAST database
  - Taxonomy file: `protein_db_taxonomy.tsv`
- **Output Files**:
  - BLAST search results
  - Taxonomy assignment results

### Assigning Taxonomy Using Genome Detective
To set a ground truth, we used Genome Detective Dengue Virus Typing Tool (https://www.genomedetective.com/app/typingtool/dengue/), which is designed to use BLAST and phylogenetic methods to identify the Dengue virus serotypes, genotypes, and major lineages of a nucleotide sequence.

- **Tools**: Genome Detective
- **Input Files**:
  - Query FASTA file: `passed_sequences.fasta`
  - BLAST database
  - Taxonomy file: `protein_db_taxonomy.tsv`
- **Output Files**:
  - BLAST search results
  - Taxonomy assignment results

### Generating Confusion Matrix: `1.8.1.whole_genome_confusion_matrix.sh`
To evaluate the accuracy of our taxonomic classification, we generated a confusion matrix using the script `1.8.1.whole_genome_confusion_matrix.sh`. This script utilized a Python script (`whole_genome_confusion_matrix.py`) to process the classification results and generate the matrix. We compared the taxonomy assigned by Genome Detective and the custom BLAST database.

- **Script**: `1.8.1.whole_genome_confusion_matrix.sh`
- **Packages and Libraries**: Python (version 3.9.15)
- **Input File**: `whole_genome_confusion.csv`
- **Output File**: Confusion matrix results

### Fetching and Indexing Dengue Virus Reference Genome – Process 2

### Fetching Reference Genomes: `2.1.fetch_ref_genome.sh`
The script `2.1.fetch_ref_genome.sh` was used to retrieve Dengue virus reference genomes from the National Center for Biotechnology Information (NCBI) database. This script filtered genomes by completeness and relevance, saving them in the specified directory for further processing.

- **Script**: `2.1.fetch_ref_genome.sh`
- **Packages and Libraries**: Entrez Direct (version 13.2)
- **Input Files**: NCBI nucleotide database query results
- **Output Files**: Dengue virus reference genomes in FASTA format

### Indexing Reference Genomes: `2.2.index_ref_genome.sh`
After fetching, the reference genomes were indexed using the script `2.2.index_ref_genome.sh`. This script utilized the Burrows-Wheeler Aligner (BWA) tool (version 0.7) to create indices for efficient sequence alignment.

- **Script**: `2.2.index_ref_genome.sh`
- **Packages and Libraries**: BWA (version 0.7)
- **Input Files**: Dengue virus reference genomes in FASTA format
- **Output Files**: Indexed reference genomes

### Dengue Virus-Specific Primers and Amplicons – Process 3

### Literature Mining for Dengue Virus-Specific Primers: `3.1.scrappaper.py`
Forked from 'ScrapPaper' by M. R. Rafsanjani. The script has been modified as follows:
- Removed the Google Scholar function.
- Added a function to scrape from both NCBI Pubmed Summary and Abstract formats.
- Added a function to merge the results.
- Added a function to extract the year from references.
- Modified the pandas write section to append new data instead of replacing it.

### Cleanning PCR Primers Sequence from the NCBI Data Result
This requires manual work. Go through the results from the '3.1.scrappaper.py' section and look for the Dengue PCR primers. Follow these guidelines:
- Only Dengue PCR Primers will be taken.
- Remove publications with no primers or where the result is N/A.
- Rename primers with no name to the first letter of every word in the publication name followed by a number.
- Remove artificial tags, probes, and RT elements.
- Replace spaces with underscores for primer names.
- Remove primer with no known genotype

- **Script**: `3.1.scrappaper.py`
- **Packages and Libraries**: Python (version 3.9.15), requests, BeautifulSoup
- **Input Files**: PubMed articles
- **Output Files**: Dengue virus primer metadata CSV file

### Converting Cleaned Primer Data to FASTA Format: `3.2.csvtofasta.sh`
The script `3.2.csvtofasta.sh` converted the cleaned CSV file containing primer data into FASTA format, preparing it for further processing.

- **Script**: `3.2.csvtofasta.sh`
- **Packages and Libraries**: Python (version 3.9.15)
- **Input Files**: Primer metadata CSV file
- **Output Files**: Primer sequences in FASTA format

### Checking Primer Positions: `3.3.check_primer_position.sh`
Forked from 'Benchmarking_PCR_ONT' by Amy Fitzpatrick. The script has been modified as follows:
- Added several error checks to ensure the accuracy and integrity of the primer positions.

- **Script**: `3.3.check_primer_position.sh`
- **Packages and Libraries**: BWA (version 0.7), samtools (version 1.10)
- **Input Files**: Primer FASTA file, reference genomes
- **Output Files**: Primer mapping positions TSV file

### Determining Primer Melting Temperatures: `3.4.fetch_temp.sh`
The `3.4.fetch_temp.sh` script calculated the minimum and maximum melting temperatures (Tm) for each primer sequence using Wallace's rule and the nearest neighbor method from an Excel file containing primer sequences. The results were then saved in a TSV file.

- **Script**: `3.4.fetch_temp.sh`
- **Packages and Libraries**: Anaconda (version 3.5.2), Biopython, pandas, openpyxl
- **Input Files**: Primer sequences
- **Output Files**: Primer Tm values

### Generating Primer Pairs

#### Conventional PCR: `3.5.1.conventional.sh`
To generate primer pairs for conventional PCR, the process involved combining forward and reverse primers based on their mapping positions, melting temperatures (Tm), and other sequence characteristics. The `3.5.1.conventional.sh` was executed to run the `conventional.py` script.

The `conventional.py` script started by calculating the GC content of each primer sequence, a crucial factor affecting primer stability and annealing temperature. The script ensured that each primer met specific conditions: length between 18 to 25 nucleotides, GC content ranging from 40% to 60%, and no homopolymer runs exceeding four nucleotides.

The script loaded mapping data from a TSV file containing primer positions and metadata, as well as additional primer metadata from `primer_metadata.tsv`, which included GC content calculations. After merging these datasets to create a comprehensive dataset, forward and reverse primers were separated based on their names ending with _F and _R.

The script iterated through the forward primers to find matching reverse primers that satisfied all predefined conditions. It ensured the primer pairs amplified regions within specified regions of interest (NS1, NS3, and NS5 of the Dengue virus genome) and that the Tm differences between forward and reverse primers were within 3°C for both minimum and maximum Tm. Additionally, the GC content difference between primers did not exceed 10%, and the length of the PCR product (amplicon) was between 300 to 500 base pairs.

Valid primer combinations were saved to `conventional_primer_combinations.tsv`, resulting in a total of 52 conventional primer combinations.

- **Script**: `3.5.1.conventional.sh`
- **Packages and Libraries**: Python (version 3.9.15)
- **Input Files**: Primer mapping positions TSV file
- **Output Files**: Primer combinations TSV file

### Semi-Nested PCR
xx
### Nested PCR
xx


### In-Silico PCR

#### Conventional PCR: `3.6.1.conv_insilico_PCR.sh`
The `3.6.1.conv_insilico_PCR.sh` script simulated conventional PCR reactions in silico using the primer pairs generated. The script used the in-silico PCR tool from Egon Ozer's GitHub repository (https://github.com/egonozer/in_silico_pcr) and a Python script to filter the resulting amplicons.

- **Script**: `3.6.1.conv_insilico_PCR.sh`
- **Packages and Libraries**: Perl (version 5.30), Python (version 3.9.15)
- **Tools**: In-silico PCR tool (https://github.com/egonozer/in_silico_pcr), custom Python script `filter_amplicons.py`
- **Input Files**: Primer combinations TSV file, Dengue virus sequences FASTA file
- **Output Files**: Filtered in-silico PCR amplicons

### Semi-Nested PCR
xx
### Nested PCR
xx

### Phylogenetic Tree Creation: `phylogenetic_tree.ipynb`
- Download nucleotide sequences from NCBI for a list of the following accession numbers:
  - DQ235145, AY323490, AF331718, AF253419, Y07863, D12937, X86784, DQ235152, DQ235151, DQ235153, AY193805, L06436, AF311056, DQ235149, U27495, X07755, L40361, DQ235144, DQ235150, KF815939, DQ235146, AY632536, AF013366, AF013375, AF013390, U88536, U87411, M93130, AF326573, KF917536, M18370, AF013384, AF161266, DQ525916, AY453411, D00246, M12294, AF013413, AY632541, AF013407, AY632545, AY632539, AF013397, AF013377, JX236040, JF895923, AY632535, DQ837642, EU707555, X03700, AY632540, DQ859056
- **Manual Work**:
  - Replace spaces with underscores in result.
  - Remove empty lines.
  - Use [Clustal Omega](https://www.ebi.ac.uk/jdispatcher/msa/muscle/) for alignment.
  - Create the phylogenetic tree using iTOL.
- **Output Files:**: 
  - `orthoflavivirus_sequences.fasta`
  - `tree.......... to be updated`


## Pre-requisite/Dependencies
The following modules, environments, and packages are needed in this project:

- **Modules**:
  - Anaconda
  - BWA (version 0.7)
  - samtools (version 1.10)

- **Environment**:
  - Conda Environment

- **Packages**:
  - Entrez Direct (EDirect, version 13.2)
  - Python libraries:
    - pandas (version 1.3.3)
    - BeautifulSoup (version 4.9.3)
    - requests (version 2.26.0)
    - Biopython (version 1.78)
    - tqdm (version 4.62.2)
    - logging (version 0.5.1.2)
    - openpyxl (version 3.0.9)
  - Perl (version 5.30)
  - VADR (version 1.5.1)
  - MMseqs2 (version 13-45111)
  - seqkit (version 0.15.0)
  - BLAST (version 2.12.0)
  - FROGS library (https://github.com/geraldinepascal/FROGS-wrappers)
  - Assign-Taxonomy-with-BLAST (https://github.com/Joseph7e/Assign-Taxonomy-with-BLAST)
  - in-silico-pcr (https://github.com/egonozer/in_silico_pcr)