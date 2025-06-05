#!/usr/bin/env bash
set -euo pipefail

########################################################################
# (1) Get the environment ready: install tools, configure EMBOSS, etc.
########################################################################

# Force the script to run relative to its own location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

#---- 1.1 Install (or confirm) required packages ----
#      emboss (includes getorf, patmatmotifs, seqret, prosextract, etc.)
#      emboss-data (for built‐in DBs, though we'll override PROSITE)
#      wget       (to fetch PROSITE files)
echo "### Installing required packages (if not already present)..."
sudo apt-get update -qq
sudo apt-get install -y emboss emboss-data wget

#---- 1.2 Set up EMBOSS_DATA for PROSITE ----
#      We’ll create a local EMBOSS_DATA directory to hold PROSITE files.
export EMBOSS_DATA="$HOME/emboss_data"
PROSITE_DIR="$EMBOSS_DATA/PROSITE"
mkdir -p "$PROSITE_DIR"

#---- 1.3 Download PROSITE if not already there ----
cd "$PROSITE_DIR"
if [[ ! -f prosite.dat || ! -f prosite.doc ]]; then
  echo "### Downloading PROSITE database..."
  wget -q https://ftp.expasy.org/databases/prosite/prosite.dat
  wget -q https://ftp.expasy.org/databases/prosite/prosite.doc
else
  echo "### PROSITE files already present, skipping download."
fi

#---- 1.4 Run prosextract so that patmatmotifs can use PROSITE ----
echo "### Running prosextract to process PROSITE..."
# prosextract will write files like prosite.pat, prosite.lines, etc., under $PROSITE_DIR
prosextract >/dev/null
cd "$SCRIPT_DIR"

########################################################################
# (2) Initialize global variables (paths for input/output)
########################################################################

# 2.1 Input GenBank file (containing one or more records; we’ll extract full sequences)
GENBANK_IN="../../data/virus.gb"

# 2.2 Intermediate/Output FASTA for whole sequence
FULLFASTA="../../data/full_sequence.fasta"

# 2.3 Intermediate FASTA of translated ORFs
ORF_FASTA="../../outputs/emboss/orfs.fasta"

# 2.4 Directory where we’ll split ORFs into individual FASTA files
ORF_SPLIT_DIR="../../outputs/emboss/orf_splits"
mkdir -p "$ORF_SPLIT_DIR"

# 2.5 Directory where per-ORF domain reports will go
DOMAIN_OUT_DIR="../../outputs/emboss/per_orf_analysis"
mkdir -p "$DOMAIN_OUT_DIR"
# 2.6 File to collect all domain‐scan results in one place
COMBINED_RESULTS="../../outputs/emboss/full_domain_analysis.txt"

# 2.7 Configuration parameters
MIN_ORF_SIZE_NT=150      # Minimum size for ORFs in nucleotides
MIN_MOTIF_LENGTH=10      # Minimum length to consider a motif "long"
MIN_HIT_COUNT=0          # Minimum number of hits required to report an ORF

if [[ ! -f "$GENBANK_IN" ]]; then
  echo "Reference file not found: $GENBANK_IN"
  exit 1
fi


########################################################################
# (3) Generate a FASTA with the entire nucleotide sequence from GenBank
########################################################################

echo "### Extracting full nucleotide sequence from GenBank..."
seqret -sequence "$GENBANK_IN" -out "$FULLFASTA" -osformat fasta

echo "   → Full-sequence FASTA written to: $FULLFASTA"

########################################################################
# (4) Get ORFs from the input FASTA
########################################################################

echo "### Extracting ORFs from $FULLFASTA..."
# getorf finds all ORFs and writes them as translated protein sequences.
getorf \
  -sequence "$FULLFASTA" \
  -outseq "$ORF_FASTA" \
  -minsize "$MIN_ORF_SIZE_NT" \
  >/dev/null

echo "   → Translated ORFs written to: $ORF_FASTA"

########################################################################
# (5) Run patmatmotifs on each extracted ORF
########################################################################

# Step 5.1: Split the multi-FASTA ORF file into individual FASTAs
echo "### Splitting $ORF_FASTA into single-ORF FASTA files..."
seqretsplit \
  -sequence "$ORF_FASTA" \
  -outseq "$ORF_SPLIT_DIR" \
  -osformat fasta \

mv ./*.fasta "$ORF_SPLIT_DIR"
echo "   → Individual ORF FASTAs in: $ORF_SPLIT_DIR"

tree $ORF_SPLIT_DIR
# Step 5.2: For each ORF FASTA, run patmatmotifs and save to DOMAIN_OUT_DIR
echo "### Running patmatmotifs on each ORF..."
# Clear combined results file (we’ll append as we go)
echo "" > "$COMBINED_RESULTS"

for orffile in "$ORF_SPLIT_DIR"/*.fasta; do
  base="$(basename "$orffile" .fasta)"
  outfile="$DOMAIN_OUT_DIR/${base}_domains.txt"

  # Run patmatmotifs: no need to specify -datafile, since EMBOSS_DATA is set
  patmatmotifs \
    -sequence "$orffile" \
    -outfile "$outfile" \
    >/dev/null

  # Extract the HitCount value from the output file
  hit_count=$(grep -m 1 "# HitCount:" "$outfile" | awk '{print $3}')
  has_long_motif=$(awk '/Length = / {if ($3 > '"$MIN_MOTIF_LENGTH"') print "yes"; exit}' "$outfile")

  # Only append to the combined results file if hitcount > minimal threshold and has long motif
  if [[ "$hit_count" -gt "$MIN_HIT_COUNT" ]] && [[ "$has_long_motif" == "yes" ]]; then
    {
      echo "### ORF: $base"
      cat "$outfile"
      echo -e "\n\n"
    } >> "$COMBINED_RESULTS"
  fi
done

echo "   → Per-ORF domain reports in: $DOMAIN_OUT_DIR"
echo "   → All domain results merged in: $COMBINED_RESULTS"

########################################################################
# (6) Done
########################################################################

echo "### Pipeline complete!"
echo "  1. Full-sequence FASTA:      $FULLFASTA"
echo "  2. Translated ORFs:          $ORF_FASTA"
echo "  3. Split ORFs:               $ORF_SPLIT_DIR"
echo "  4. Per-ORF domain reports:   $DOMAIN_OUT_DIR"
echo "  5. Combined domain results:  $COMBINED_RESULTS"
