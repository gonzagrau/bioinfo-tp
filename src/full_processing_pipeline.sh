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

# TODO: Prepare the environment and run stuff from other scripts at 01, 02,