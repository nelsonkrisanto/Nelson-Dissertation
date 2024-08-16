# Development of an In-Silico Pipeline for Assessing Specificity of Dengue Virus Genotyping via Conventional PCR Assay

## Introduction

This repository contains the resources, scripts, and data related to the dissertation project titled **"Development of an In-Silico Pipeline for Assessing Specificity of Dengue Virus Genotyping via Conventional PCR Assay"** by Nelson Krisanto.

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
