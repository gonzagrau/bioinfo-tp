import os
from Bio import SeqIO, AlignIO
import subprocess

# Rutas a los ejecutables (ajusta según tu sistema)
muscle_exe = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/Lou_campeona/muscle-win64.v5.3.exe"
mafft_exe = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/Lou_campeona/mafft-win/mafft.bat"

fastas_folder = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/output/fastas"
output_muscle = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/Lou_campeona/output_muscle"
output_mafft = r"C:/Users/lourd/OneDrive/Documentos/00_bioinfo/bioinfo-tp/Lou_campeona/output_mafft"
os.makedirs(output_muscle, exist_ok=True)
os.makedirs(output_mafft, exist_ok=True)

tmp_prot = os.path.join(output_muscle, "all_proteins.fasta")
tmp_dna = os.path.join(output_muscle, "all_dna.fasta")
tmp_prot_mafft = os.path.join(output_mafft, "all_proteins.fasta")
tmp_dna_mafft = os.path.join(output_mafft, "all_dna.fasta")

def is_dna(seq, threshold=0.95):
    """Devuelve True si la secuencia parece ADN."""
    seq = str(seq).upper()
    dna_letters = set("ACGTN")
    if len(seq) == 0:
        return False
    return sum(1 for c in seq if c in dna_letters) / len(seq) >= threshold

# 1. Clasifica y guarda las secuencias para MUSCLE y MAFFT
with open(tmp_prot, "w") as prot_out, open(tmp_dna, "w") as dna_out, \
     open(tmp_prot_mafft, "w") as prot_out_mafft, open(tmp_dna_mafft, "w") as dna_out_mafft:
    for fname in os.listdir(fastas_folder):
        if fname.endswith(".fasta") or fname.endswith(".fa"):
            for record in SeqIO.parse(os.path.join(fastas_folder, fname), "fasta"):
                record_id = f"{os.path.splitext(fname)[0]}|{record.id}"
                record.description = ""
                record.id = record_id
                if is_dna(record.seq):
                    SeqIO.write(record, dna_out, "fasta")
                    SeqIO.write(record, dna_out_mafft, "fasta")
                else:
                    SeqIO.write(record, prot_out, "fasta")
                    SeqIO.write(record, prot_out_mafft, "fasta")

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

def run_mafft(input_fasta, output_fasta):
    if os.path.exists(input_fasta) and os.path.getsize(input_fasta) > 0:
        cmd = [mafft_exe, "--auto", input_fasta]
        print("Ejecutando:", " ".join(cmd))
        with open(output_fasta, "w") as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print(f"Error ejecutando MAFFT para {input_fasta}:")
            print(result.stderr)
            return None
        return output_fasta
    else:
        print(f"No hay secuencias en {input_fasta}, no se ejecuta MAFFT.")
        return None

# 2. Ejecuta MUSCLE y MAFFT para proteínas y ADN
out_prot_muscle = os.path.join(output_muscle, "aln_proteins.fasta")
out_dna_muscle = os.path.join(output_muscle, "aln_dna.fasta")
out_prot_mafft = os.path.join(output_mafft, "aln_proteins.fasta")
out_dna_mafft = os.path.join(output_mafft, "aln_dna.fasta")

prot_aln_muscle = run_muscle(tmp_prot, out_prot_muscle)
dna_aln_muscle = run_muscle(tmp_dna, out_dna_muscle)
prot_aln_mafft = run_mafft(tmp_prot_mafft, out_prot_mafft)
dna_aln_mafft = run_mafft(tmp_dna_mafft, out_dna_mafft)

# 3. Imprime los alineamientos
for metodo, prot_aln, dna_aln in [
    ("MUSCLE", prot_aln_muscle, dna_aln_muscle),
    ("MAFFT", prot_aln_mafft, dna_aln_mafft)
]:
    for label, aln_file in [("PROTEÍNAS", prot_aln), ("ADN", dna_aln)]:
        if aln_file and os.path.exists(aln_file):
            print(f"\nAlineamiento {metodo} ({label}):")
            alignment = AlignIO.read(aln_file, "fasta")
            print(f"Número de secuencias alineadas: {len(alignment)}")
            print(f"Longitud del alineamiento: {alignment.get_alignment_length()}")
            for record in alignment:
                print(f">{record.id}\n{record.seq}\n")