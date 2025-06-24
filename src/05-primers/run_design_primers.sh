#!/usr/bin/env bash
set -euo pipefail

# Resolve paths relative to this script's location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTA_DIR="$SCRIPT_DIR/../../outputs/fastas"
PARAMS="$SCRIPT_DIR/params.json"
OUT_DIR="$SCRIPT_DIR/../../outputs/primers"

mkdir -p "$OUT_DIR"

for fasta in "$FASTA_DIR"/*_cds.fasta; do
    echo "Diseñando primers para $(basename "$fasta")"
    python "$SCRIPT_DIR/design_primers.py" "$fasta" "$PARAMS" "$OUT_DIR"
done

echo "Resultados del diseño de primers guardados en $OUT_DIR"


