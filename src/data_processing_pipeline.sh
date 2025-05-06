#!/bin/bash

# get environment ready
pip install -r ../requirements.txt
if [[ ! -d ../data/swissprot/ ]]; then
	mkdir ../data/swissprot/
	cd ../data/swissprot/
	update_blastdb.pl --decompress swissprot
	cd ../src/
	
# get fastas
python ./read_seq.py

# run blast
blastp -query ../output/fastas/lep-sequence_protein_sequences.fasta \
       -db ../data/swissprot/swissprot \
       -out ../output/fasta_results.txt \
       -outfmt 6
