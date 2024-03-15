# Title
## High Level Guide to Using the Scripts

### Configuration: `00.config.sh`
Please update the script with your specific settings where you see the comment "## adjust". Ensure these settings align with your project requirements.

### Fetching NCBI Data: `01.fetch_NCBI.sh`
Forked from 'Benchmarking_PCR_ONT' by Amy Fitzpatrick. The script has been modified as follows:
Modify the sections marked as "## adjust" to suit your research needs. When you reach the portion of the script that says, "# Perform the search to get dengue virus sequences. Adjust the search term as needed," refer to the instructions provided below.

To tailor the search query for your need, use the NCBI Nucleotide Advanced Search Builder. The image below (Figure 1.1) illustrates the Advanced Search Builder interface, which you can access at the following URL: [NCBI Nucleotide Advanced Search Builder](https://www.ncbi.nlm.nih.gov/nuccore/advanced).

![Figure 1.1 - NCBI Nucleotide Advanced Search Builder](images/01_01.png)
<p align="center">
  Figure 1.1 - NCBI Nucleotide Advanced Search Builder
</p>

In the Advanced Search Builder, adjust the search parameters as necessary. Your configured parameters will reflect in the search details field, as highlighted in Figure 1.1.

### Fetching NCBI Data: `02.scrappaper.py`
Forked from 'ScrapPaper' by M. R. Rafsanjani. The script has been modified as follows:
- Removed the Google Scholar function.
- Added a function to scrape from both NCBI Pubmed Summary and Abstract formats.
- Added a function to merge the results.
- Added a function to extract the year from references.
- Modified the pandas write section to append new data instead of replacing it.

### Scrape PCR Primers Sequence from the NCBI Data Result
This requires manual work. Go through the results from the 'Fetching NCBI Data' section and look for the Dengue PCR primers. Follow these guidelines:
- Only Dengue PCR Primers will be taken.
- Remove publications with no primers or where the result is N/A.
- Rename primers with no name to the first letter of every word in the publication name followed by a number.
- Remove artificial tags, probes, and RT elements.
- Replace spaces with underscores for primer names.

### Convert cleanned CSV PCR Primers result to FASTA format: `03.csvtofasta.py`
.......to be updated

### Check Primer Position: `04.check_primer_position.sh`
Forked from 'Benchmarking_PCR_ONT' by Amy Fitzpatrick. The script has been modified as follows:
- Added several error check
- 
.......to be updated



## Pre-requisite/Dependencies
The following module, environment, and packages are needed in this project:
* Module: 
  * Anaconda
* Environment:
  * Conda Environment
* Packages:
  * Entrez Direct (EDirect)

tes