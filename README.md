# KUZM1124

![version](https://img.shields.io/badge/current_v-1.0.0-blue)
![version](https://img.shields.io/badge/python-3.8+-blue)

## Description
KUZM1124 is a bioinformatics pipeline designed to process raw Illumina sequencing reads and verify them against a library database. 
It utilizes various tools and scripts to ensure the integrity and accuracy of the sequencing data.

## Features
- **Input**: Accepts raw Illumina sequencing reads in FASTQ format.
- **Merge Reads**: Merges paired-end reads using `PandaSeq`.
- **Quality Control**: Performs quality control checks using `FastQC` and `MultiQC`.
- **Database Verification**: Compares processed reads against a library database by direct search.
- **Output**: Generates a report summarizing the results of the verification process.

## Requirements
- Python 3.8 or higher
- PandaSeq
- FastQC
- MultiQC
- Other dependencies as specified in the `requirements.txt` file.

## Installation
- Create a virtual environment:
  ```bash
  python -m venv kuz_env
  source kuz_env/bin/activate  # On Windows use: kuz_env\Scripts\activate
  ```
- Install the required packages:
  ```bash
    pip install -r requirements.txt
  ```