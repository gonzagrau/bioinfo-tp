#!/bin/bash
set -e # makes the script exit immediately if any command returns a non-zero 

# 1. Prepare environment
pip install -r ../requirements.txt

# 2. Download Swiss‐Prot DB if missing
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


# 5. Per‐sequence BLASTp for NP_000221.1
blastp \
  -query ../output/fastas/lep-sequence_NP_000221.1_protein.fasta \
  -db    ../data/swissprot/swissprot \
  -out   ../output/results_NP_000221.1.tsv \
  -evalue 1e-5 \
  -outfmt '6 qseqid sseqid pident length evalue bitscore'

#    prepend header + align columns
sed '1i qseqid\tsseqid\tpident\tlength\tevalue\tbitscore' \
  ../output/results_NP_000221.1.tsv | \
  column -t -s $'\t' > ../output/results_NP_000221.1.tmp.tsv
mv ../output/results_NP_000221.1.tmp.tsv ../output/results_NP_000221.1.tsv

# 6. Per‐sequence BLASTp for XP_005250397.1
blastp \
  -query ../output/fastas/lep-sequence_XP_005250397.1_protein.fasta \
  -db    ../data/swissprot/swissprot \
  -out   ../output/results_XP_005250397.1.tsv \
  -evalue 1e-5 \
  -outfmt '6 qseqid sseqid pident length evalue bitscore'

#    prepend header + align columns
sed '1i qseqid\tsseqid\tpident\tlength\tevalue\tebit sscore' \
  ../output/results_XP_005250397.1.tsv | \
  column -t -s $'\t' > ../output/results_XP_005250397.1.tmp.tsv
mv ../output/results_XP_005250397.1.tmp.tsv ../output/results_XP_005250397.1.tsv
