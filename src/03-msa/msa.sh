#!/bin/bash

# MSA Bash Script
# This script runs multiple sequence alignment tools (MUSCLE, MAFFT, Clustal Omega)
# on DNA and protein sequences.
# Ensure the script is run from the correct directory

cd "$(dirname "$0")" || exit 1
# Install required tools if not already installed
if ! command -v muscle &> /dev/null; then
    echo "MUSCLE not found, installing..."
    sudo apt-get install muscle -y
else
    echo "MUSCLE is already installed."
fi
if ! command -v mafft &> /dev/null; then
    echo "MAFFT not found, installing..."
    sudo apt-get install mafft -y
else
    echo "MAFFT is already installed."
fi
if ! command -v clustalo &> /dev/null; then
    echo "Clustal Omega not found, installing..."
    sudo apt-get install clustalo -y
else
    echo "Clustal Omega is already installed."
fi

# Set up logging
echo "Starting MSA script..."

# Define input and output directories
BASE_OUTPUT_DIR="../../outputs/multiple-sequence-alignment"
mkdir -p "$BASE_OUTPUT_DIR"
OUTPUT_MUSCLE="$BASE_OUTPUT_DIR/output_muscle"
OUTPUT_MAFFT="$BASE_OUTPUT_DIR/output_mafft"
OUTPUT_CLUSTALO="$BASE_OUTPUT_DIR/output_clustalo"

# Create output directories
echo "Creating output directories..."
mkdir -p "$OUTPUT_MUSCLE"
mkdir -p "$OUTPUT_MAFFT"
mkdir -p "$OUTPUT_CLUSTALO"
echo "Directories created."

# Input files
INPUT_DNA="../../outputs/fastas/joined/all_cds.fasta"
INPUT_PROT="../../outputs/fastas/joined/all_proteins.fasta"

# Output files
OUT_PROT_MUSCLE="$OUTPUT_MUSCLE/aln_proteins.fasta"
OUT_DNA_MUSCLE="$OUTPUT_MUSCLE/aln_dna.fasta"
OUT_PROT_MAFFT="$OUTPUT_MAFFT/aln_proteins.fasta"
OUT_DNA_MAFFT="$OUTPUT_MAFFT/aln_dna.fasta"
OUT_PROT_CLUSTALO="$OUTPUT_CLUSTALO/aln_proteins.fasta"
OUT_DNA_CLUSTALO="$OUTPUT_CLUSTALO/aln_dna.fasta"

# Run MUSCLE
run_muscle() {
    local input=$1
    local output=$2
    echo "Starting MUSCLE alignment for $(basename $input)"
    
    # Check if input file exists and is not empty
    if [ -f "$input" ] && [ -s "$input" ]; then
        muscle -in "$input" -out "$output"
        if [ $? -eq 0 ]; then
            echo "MUSCLE executed successfully, results in: $output"
        else
            echo "Error executing MUSCLE for $input"
        fi
    else
        echo "Input file $input not found or empty, skipping MUSCLE alignment"
    fi
}

# Run MAFFT
run_mafft() {
    local input=$1
    local output=$2
    echo "Starting MAFFT alignment for $(basename $input)"
    
    # Check if input file exists and is not empty
    if [ -f "$input" ] && [ -s "$input" ]; then
        mafft --auto "$input" > "$output"
        
        if [ $? -eq 0 ]; then
            echo "MAFFT executed successfully, results in: $output"
        else
            echo "Error executing MAFFT for $input"
        fi
    else
        echo "Input file $input not found or empty, skipping MAFFT alignment"
    fi
}

# Run Clustal Omega
run_clustalo() {
    local input=$1
    local output=$2
    echo "Starting Clustal Omega alignment for $(basename $input)"
    
    # Check if input file exists and is not empty
    if [ -f "$input" ] && [ -s "$input" ]; then
        clustalo -i "$input" -o "$output" --force --outfmt=fasta
        
        if [ $? -eq 0 ]; then
            echo "Clustal Omega executed successfully, results in: $output"
        else
            echo "Error executing Clustal Omega for $input"
        fi
    else
        echo "Input file $input not found or empty, skipping Clustal Omega alignment"
    fi
}

# Execute alignments
echo "Running alignments..."

# Run MUSCLE alignments
run_muscle "$INPUT_PROT" "$OUT_PROT_MUSCLE"
run_muscle "$INPUT_DNA" "$OUT_DNA_MUSCLE"

# Run MAFFT alignments
run_mafft "$INPUT_PROT" "$OUT_PROT_MAFFT"
run_mafft "$INPUT_DNA" "$OUT_DNA_MAFFT"

# Run Clustal Omega alignments
run_clustalo "$INPUT_PROT" "$OUT_PROT_CLUSTALO"
run_clustalo "$INPUT_DNA" "$OUT_DNA_CLUSTALO"

echo "All alignments completed."