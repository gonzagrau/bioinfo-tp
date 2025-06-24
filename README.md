# Bioinfo-TP

Final project for our Bioinformatics class

## Overview
This project provides tools for processing GenBank (.gb) DNA sequence files. It extracts coding sequences (CDS), translates them, and stores them in various formats. It then runs BLASTp on the protein sequences, and applies MSA to the relevant results. There is also an EMBOSS pipeline to extract motifs from ORFs, starting back from the raw .gb, and finally a primer design module.

## Requirements
- Python 3.6+

Install requirements with:
```bash
pip install -r requirements.txt
```

## Features at `src`
- `src/01-preprocess`: Extract ORFs from genbank file and translate them
- `src/02-BLAST`: Run BLASTP on them
- `src/03-msa`: Run multpiple sequence alignment on those
- `src/04-emboss`: Find motifs on ORFs using EMBOSS
- `src/05-primers`: Design primers for the given sequence

## Usage
To run the entire analysis (preprocessing, BLAST search and multiple sequence
alignment) in one command, execute the `run_pipeline.sh` script from the project
root:

```bash
bash full_pipeline_script.sh
```

The script installs the required Python packages and sequentially runs:

1. `src/01-preprocess/read_seq.py`
2. `src/02-BLAST/blast_script.sh`
3. `src/03-msa/msa.sh`

Results will be stored inside the `outputs` directory.

To run the `src/05-primers` use the `src/05-primers/run_design_primers.sh` bash script, 
it installs the required libraries for the Python file and also iterates over the Fasta sequences.
