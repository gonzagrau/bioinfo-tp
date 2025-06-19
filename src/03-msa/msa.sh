#!/bin/bash

# MSA Bash Script
# This script runs multiple sequence alignment tools (MUSCLE, MAFFT, Clustal Omega)
# on DNA and protein sequences.

# Set up logging
echo "Starting MSA script..."

# Find binaries
find_binary() {
    local name=$1
    local path=$(which $name)
    if [ -n "$path" ]; then
        echo $path
        return
    fi
    
    # Common Linux binary locations to check
    local locations=(
        "/usr/bin/$name"
        "/usr/local/bin/$name"
        "/opt/$name/$name"
    )
    
    for loc in "${locations[@]}"; do
        if [ -f "$loc" ] && [ -x "$loc" ]; then
            echo $loc
            return
        fi
    done
    
    # Use the name and hope it's in PATH
    echo $name
}

echo "Looking for alignment binaries..."
MUSCLE_BIN=$(find_binary "muscle")
MAFFT_BIN=$(find_binary "mafft")
CLUSTALO_BIN=$(find_binary "clustalo")
echo "Found binaries: MUSCLE: $MUSCLE_BIN, MAFFT: $MAFFT_BIN, Clustal Omega: $CLUSTALO_BIN"

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
        # Check MUSCLE version to determine command format
        if [[ "$MUSCLE_BIN" == *"5"* ]] || [[ "$($MUSCLE_BIN -h 2>&1)" == *"-align"* ]]; then
            $MUSCLE_BIN -align "$input" -output "$output"
        else
            # Older MUSCLE versions
            $MUSCLE_BIN -in "$input" -out "$output"
        fi
        
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
        $MAFFT_BIN --auto "$input" > "$output"
        
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
        $CLUSTALO_BIN -i "$input" -o "$output" --force --outfmt=fasta
        
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