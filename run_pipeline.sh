#!/bin/bash
set -e

# This script installs dependencies and runs the processing
# steps located in src/01-preprocess, src/02-BLAST and src/03-msa.

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

# 1. Install Python requirements
pip install -r "$ROOT_DIR/requirements.txt"

# 2. Run preprocessing (extract CDS and FASTA files)
cd "$ROOT_DIR/src/01-preprocess"
python read_seq.py

# 3. Run BLAST processing on generated FASTA files
cd "$ROOT_DIR/src/02-BLAST"
bash blast_script.sh

# 4. Run multiple sequence alignment
cd "$ROOT_DIR/src/03-msa"
bash msa.sh

# Return to root directory at the end
cd "$ROOT_DIR"
