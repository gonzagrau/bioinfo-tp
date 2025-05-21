import os
from Bio import SeqIO, AlignIO
import subprocess

# Rutas a los ejecutables (ajusta según tu sistema)
muscle_exe = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/Lou_campeona/muscle-win64.v5.3.exe"

fastas_folder = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/output/fastas"
output_folder = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/Lou_Campeona/output"
os.makedirs(output_folder, exist_ok=True)

tmp_prot = os.path.join(output_folder, "all_proteins.fasta")
tmp_dna = os.path.join(output_folder, "all_dna.fasta")

def is_dna(seq, threshold=0.95):
    """Devuelve True si la secuencia parece ADN."""
    seq = str(seq).upper()
    dna_letters = set("ACGTN")
    if len(seq) == 0:
        return False
    return sum(1 for c in seq if c in dna_letters) / len(seq) >= threshold

# 1. Clasifica y guarda las secuencias
with open(tmp_prot, "w") as prot_out, open(tmp_dna, "w") as dna_out:
    for fname in os.listdir(fastas_folder):
        if fname.endswith(".fasta") or fname.endswith(".fa"):
            for record in SeqIO.parse(os.path.join(fastas_folder, fname), "fasta"):
                record.id = f"{os.path.splitext(fname)[0]}|{record.id}"
                record.description = ""
                if is_dna(record.seq):
                    SeqIO.write(record, dna_out, "fasta")
                else:
                    SeqIO.write(record, prot_out, "fasta")

def run_muscle(input_fasta, output_fasta):
    if os.path.exists(input_fasta) and os.path.getsize(input_fasta) > 0:
        cmd = [muscle_exe, "-align", input_fasta, "-output", output_fasta]
        print("Ejecutando:", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error ejecutando MUSCLE para {input_fasta}:")
            print(result.stderr)
            return None
        return output_fasta
    else:
        print(f"No hay secuencias en {input_fasta}, no se ejecuta MUSCLE.")
        return None

# 2. Ejecuta MUSCLE para proteínas y ADN
out_prot = os.path.join(output_folder, "aln_proteins.fasta")
out_dna = os.path.join(output_folder, "aln_dna.fasta")

prot_aln = run_muscle(tmp_prot, out_prot)
dna_aln = run_muscle(tmp_dna, out_dna)

# 3. Imprime los alineamientos
for label, aln_file in [("PROTEÍNAS", prot_aln), ("ADN", dna_aln)]:
    if aln_file and os.path.exists(aln_file):
        print(f"\nAlineamiento MUSCLE ({label}):")
        alignment = AlignIO.read(aln_file, "fasta")
        print(f"Número de secuencias alineadas: {len(alignment)}")
        print(f"Longitud del alineamiento: {alignment.get_alignment_length()}")
        for record in alignment:
            print(f">{record.id}\n{record.seq}\n")