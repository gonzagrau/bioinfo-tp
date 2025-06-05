#!/bin/bash
set -e # makes the script exit immediately if any command returns a non-zero 

# 1. Prepare environment
pip install -r ../requirements.txt

# 2. Install BLAST+ if not installed
if ! command -v blastp &> /dev/null; then
  echo "BLAST+ not found. Installing..."
  if command -v apt-get &> /dev/null; then
    sudo apt-get update
    sudo apt-get install -y ncbi-blast+
  elif command -v yum &> /dev/null; then
    sudo yum install -y ncbi-blast+
  else
    echo "Please install BLAST+ manually."
    exit 1
  fi
fi

# 3. Download Swiss‐Prot DB if missing
if [[ ! -d ../data/swissprot ]]; then
  mkdir -p ../data/swissprot
  cd ../data/swissprot
  update_blastdb.pl --decompress swissprot
  cd - #>/dev/null
fi

# 3. Generate FASTAs
python ../src/read_seq.py

# 4. Make sure output dirs exist
mkdir -p ../output
mkdir -p ../output/fastas


# 5. Per‐sequence BLASTp for each FASTA in ../output/fastas/
for fasta in ../output/fastas/*.fasta; do
  base=$(basename "$fasta" .fasta)
  out_tsv="../output/results_${base}.tsv"
  tmp_tsv="../output/results_${base}.tmp.tsv"

  blastp \
    -query "$fasta" \
    -db ../data/swissprot/swissprot \
    -out "$out_tsv" \
    -evalue 1e-5 \
    -outfmt '6 qseqid sseqid pident length evalue bitscore'

  # prepend header + align columns
  sed '1i qseqid\tsseqid\tpident\tlength\tevalue\tbitscore' "$out_tsv" | \
    column -t -s $'\t' > "$tmp_tsv"
  mv "$tmp_tsv" "$out_tsv"
done



