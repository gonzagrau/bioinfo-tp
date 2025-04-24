# Bioinfo-TP: GenBank CDS Extractor

Final project for our Bioinformatics class

## Overview
This project provides tools for processing GenBank (.gb) DNA sequence files. It extracts coding sequences (CDS), translates them, and stores them in various formats.

## Features
- Extract CDS features from GenBank files
- Get the nucleotide sequence for each CDS
- Translate CDS sequences to protein sequences
- Store data in CSV and FASTA formats

## Usage
Run the script from the src directory:

```bash
cd src
python read_seq.py
```

This will:
1. Read the GenBank file from `../data/lep-sequence.gb`
2. Extract all CDS features
3. Print information about each CDS
4. Save the data to the `output` directory:
   - `cds_features.csv`: All CDS information in CSV format
   - `cds_sequences.fasta`: CDS nucleotide sequences in FASTA format
   - `protein_sequences.fasta`: Translated protein sequences in FASTA format

## Requirements
- Python 3.6+
- Biopython
- pandas

Install requirements with:
```bash
pip install -r requirements.txt
```
