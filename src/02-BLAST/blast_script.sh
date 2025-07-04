#!/bin/bash
set -e # makes the script exit immediately if any command returns a non-zero

# Install BLAST+ if not already present
if ! command -v blastp >/dev/null 2>&1; then
    echo "### Installing NCBI BLAST+..."
    if command -v apt-get >/dev/null 2>&1; then
        apt-get update -y
        apt-get install -y ncbi-blast+
    else
        echo "Error: apt-get not available. Please install NCBI BLAST+ manually." >&2
        exit 1
    fi
fi

# Prepare environment
SWISSPROT_DIR="../../data/swissprot"
if [[ ! -d "$SWISSPROT_DIR" ]]; then
  echo "### Downloading Swiss‐Prot database..."
  mkdir -p "$SWISSPROT_DIR"
  cd "$SWISSPROT_DIR"
  update_blastdb --decompress swissprot
  cd -
else
  echo "### Swiss‐Prot database already present."
fi
export BLASTDB="$(realpath "$SWISSPROT_DIR")"

# Make sure output dirs exist
mkdir -p ../../outputs
mkdir -p ../../outputs/fastas
mkdir -p ../../outputs/blast_results


# Process all FASTA files in the fastas directory
for fasta in ../../outputs/fastas/*_protein.fasta; do
    # Extract the accession number from the filename
    filename=$(basename "$fasta")
    accession="${filename#*_}" # Get the part after the first underscore
    # Run BLASTp

    echo "Processing sequence $fasta..."
    blastp \
            -query "$fasta" \
      -db    swissprot \
            -out   "../../outputs/blast_results/results_${accession}.tsv" \
      -evalue 1e-5 \
      -outfmt '6 qseqid sseqid pident length evalue bitscore'

        # Prepend header + align columns
    sed '1i qseqid\tsseqid\tpident\tlength\tevalue\tbitscore' \
        "../../outputs/blast_results/results_${accession}.tsv" | \
            column -t -s $'\t' > "../../outputs/blast_results/results_${accession}.tmp.tsv"
        mv "../../outputs/blast_results/results_${accession}.tmp.tsv" "../../outputs/blast_results/results_${accession}.tsv"

    echo "BLASTp results for $accession saved to ../../outputs/blast_results/results_${accession}.tsv"
done

echo "### BLASTp processing completed for all sequences."
