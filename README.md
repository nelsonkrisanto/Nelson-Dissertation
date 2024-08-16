
# Development of an In-Silico Pipeline for Assessing Specificity of Dengue Virus Genotyping via Conventional PCR Assay

## Introduction

This repository contains the resources, scripts, and data related to the dissertation project titled **"Development of an In-Silico Pipeline for Assessing Specificity of Dengue Virus Genotyping via Conventional PCR Assay"**.

The project aims to develop and evaluate an in-silico pipeline for assessing the specificity of dengue virus genotyping using conventional PCR assays. The primary objective is to establish a reliable and efficient method for identifying dengue virus serotypes and genotypes from PCR-amplified sequences, which is critical for accurate diagnosis and epidemiological studies.

## Table of Contents

- [Introduction](#introduction)
- [Project Overview](#project-overview)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Description](#pipeline-description)
- [Results](#results)
- [Figures and Tables](#figures-and-tables)
- [References](#references)
- [Contact](#contact)

## Project Overview

The dissertation focuses on the development of an in-silico pipeline that:

1. **Fetches and processes dengue virus sequences** from public databases.
2. **Clusters sequences** to reduce redundancy and improve analysis efficiency.
3. **Generates and evaluates primers** specific to different dengue virus serotypes and genotypes.
4. **Performs in-silico PCR** to simulate and assess the effectiveness of the designed primers.
5. **Analyzes the results** to determine the specificity and accuracy of the genotyping process.

## Installation

To replicate the pipeline, you will need to have the following software installed:

- **Anaconda** (Python 3.7+)
- **MMseqs2** (for sequence clustering)
- **BLAST** (for sequence alignment)
- **VADR** (for amplicon generation)
- **Biopython** (for sequence handling and analysis)
- **Pandas** (for data manipulation)
- **Matplotlib and Seaborn** (for data visualization)

Clone the repository and set up the environment:

```bash
git clone https://github.com/nelsonkrisanto/Nelson-Dissertation.git
cd Nelson-Dissertation
conda create --name dengue_env python=3.9
conda activate dengue_env
pip install -r requirements.txt
```

## Usage

The pipeline consists of several steps, each executable via the provided scripts. Below are the steps and commands:

1. **Fetch Dengue Virus Sequences:**

   This script fetches dengue virus sequences from the NCBI database.

   ```sh
   bash scripts/1.1.fetch_NCBI.sh
   ```

2. **Cluster Sequences:**

   This step clusters the downloaded sequences to reduce redundancy.

   ```sh
   bash scripts/1.2.cluster_MMseqs2.sh
   ```

3. **Fetch Orthoflavivirus Sequences and Combine:**

   Fetches sequences from the Orthoflavivirus genus and combines them with the dengue virus sequences.

   ```sh
   bash scripts/1.3.fetch_combine_orthoflavivirus.sh
   ```

4. **Generate Primers and Perform In-Silico PCR:**

   Designs primers and performs an in-silico PCR to simulate the amplification process.

   ```sh
   bash scripts/3.5.conventional.sh
   bash scripts/3.6.conv_insilico_PCR.sh
   ```

5. **Analyze Results:**

   Processes the PCR results and generates confusion matrices to evaluate the performance of the primers.

   ```sh
   python scripts/generate_confusion_matrix.py
   ```

## Pipeline Description

The pipeline is a comprehensive workflow designed to analyze dengue virus sequences and assess the specificity of genotyping using conventional PCR. It consists of the following steps:

1. **Sequence Retrieval:**

   Sequences are fetched from the NCBI database using custom scripts. These sequences are then pre-processed and filtered to ensure they meet specific criteria (e.g., sequence length).

2. **Clustering:**

   Sequences are clustered using MMseqs2 to reduce redundancy. This step is crucial for efficient downstream analysis, ensuring that only representative sequences are used.

3. **Primer Design:**

   Primers are designed for different dengue virus serotypes and genotypes. The script takes into account factors like primer length, GC content, melting temperature, and specificity.

4. **In-Silico PCR:**

   The designed primers are tested in an in-silico PCR simulation. This step helps to identify the primers that can effectively amplify the target regions of the dengue virus genome.

5. **Taxonomy Assignment:**

   Sequences are assigned taxonomy using tools like BLAST and Genome Detective. This helps in validating the specificity and accuracy of the designed primers.

6. **Results Analysis:**

   The results of the in-silico PCR are analyzed to generate confusion matrices, which help in understanding the performance of the primers in correctly identifying the serotypes and genotypes.

## Results

The results from the pipeline provide insights into the effectiveness of the designed primers:

1. **Distribution of Dengue Virus Primers:**

   The analysis shows the distribution of primers over the years and across different dengue virus serotypes. It highlights trends in primer development and the focus of research on specific serotypes.

2. **Confusion Matrices:**

   Confusion matrices are generated for both serotype and genotype classifications. These matrices illustrate how well the primers distinguish between different serotypes and genotypes, with metrics such as precision, recall, and F1 scores.

3. **Heatmaps:**

   Heatmaps for missed and mismatched classifications provide a visual representation of the accuracy of the primers. These visualizations help identify areas where primer design might need improvement.

4. **Phylogenetic Analysis:**

   The phylogenetic analysis shows the evolutionary relationships among different dengue virus strains, helping to understand the genetic diversity and its implications for primer design.

5. **Performance Metrics:**

   Detailed performance metrics, including F1 scores for each primer pair, are calculated to assess the accuracy and reliability of the PCR assays.

## Figures and Tables

Key figures and tables are included in the report to illustrate the findings, such as:

- **Figure 6:** Distribution of Dengue Virus Primers by Year and Serotype.
- **Table 4:** DENV Infectivity by Genotype.
- **Figure 13:** Serotype Confusion Matrix.
- **Figure 15:** Missed Classification Heatmap.

## References

The project references various studies and datasets, which are documented in the dissertation report. Please refer to the `References` section in the report for detailed citations.

## Contact

For any questions or collaboration requests, please contact Nelson Krisanto at [nelsonkrisanto@example.com](mailto:nelsonkrisanto@example.com).
