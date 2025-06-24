#!/bin/bash
set -e

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

#Install Python requirements
pip install -r "$ROOT_DIR/requirements.txt"

# Run preprocessing 
cd "$ROOT_DIR/src/01-preprocess"
python read_seq.py

# Run BLAST 
cd "$ROOT_DIR/src/02-BLAST"
chmod +x blast_script.sh
bash blast_script.sh

# Run multiple sequence alignment
cd "$ROOT_DIR/src/03-msa"
chmod +x msa.sh
bash msa.sh

# Return to root directory 
cd "$ROOT_DIR"
