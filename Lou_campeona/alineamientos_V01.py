import os
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MuscleCommandline
import subprocess

muscle_exe = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/Lou_campeona/muscle-win64.v5.3.exe"
fastas_folder = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/output/fastas"
tmp_fasta = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/Lou_campeona/all_seqs_tmp.fasta"
out_file = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/Lou_campeona/out_file.fasta"

# 1. Junta todas las secuencias de todos los .fasta en un solo archivo temporal
with open(tmp_fasta, "w") as outfile:
    for fname in os.listdir(fastas_folder):
        if fname.endswith(".fasta") or fname.endswith(".fa"):
            for record in SeqIO.parse(os.path.join(fastas_folder, fname), "fasta"):
                # Puedes modificar el id para que incluya el nombre del archivo (especie)
                record.id = f"{os.path.splitext(fname)[0]}|{record.id}"
                record.description = ""
                SeqIO.write(record, outfile, "fasta")

# 2. Ejecuta MUSCLE usando subprocess (más robusto para MUSCLE v5+)
cmd = [
    muscle_exe,
    "-align", tmp_fasta,
    "-output", out_file
]
print("Ejecutando:", " ".join(cmd))
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode != 0:
    print("Error ejecutando MUSCLE:")
    print(result.stderr)
    exit(1)

# 3. Lee el alineamiento resultante
alignment = AlignIO.read(out_file, "fasta")
print(f"Número de secuencias alineadas: {len(alignment)}")
print(f"Longitud del alineamiento: {alignment.get_alignment_length()}")
print("\nAlineamiento múltiple:\n")
for record in alignment:
    print(f">{record.id}\n{record.seq}\n")