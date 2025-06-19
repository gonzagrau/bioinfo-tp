import os
import subprocess
import shutil
import logging
from Bio import SeqIO, AlignIO

# For Linux-specific paths, no Windows checks needed
def find_binary(name):
    path = shutil.which(name)
    if path:
        return path
    # Common Linux binary locations to check
    locations = [
        f"/usr/bin/{name}",
        f"/usr/local/bin/{name}",
        f"/opt/{name}/{name}"
    ]
    for loc in locations:
        if os.path.isfile(loc) and os.access(loc, os.X_OK):
            return loc
    return name  # Use the name and hope it's in PATH

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Locate binaries (standard locations on Debian)
print("Buscando binarios de alineamiento...")
muscle_bin = find_binary("muscle")
mafft_bin = find_binary("mafft")
clustalo_bin = find_binary("clustalo")
print(f"Found binaries: MUSCLE: {muscle_bin}, MAFFT: {mafft_bin}, Clustal Omega: {clustalo_bin}")


def is_dna(seq, threshold=0.95):
    """Devuelve True si la secuencia parece ADN."""
    seq = str(seq).upper()
    dna_letters = set("ACGTN")
    if len(seq) == 0:
        return False
    return sum(1 for c in seq if c in dna_letters) / len(seq) >= threshold


def run_alignment(algorithm_name, binary, input_fasta, output_fasta, get_command_fn):
    """
    Generic function to run alignment tools.

    Args:
        algorithm_name: Name of the algorithm (for logging)
        binary: Path to the binary
        input_fasta: Path to input FASTA file
        output_fasta: Path for output aligned FASTA
        get_command_fn: Function that returns the command list for this algorithm

    Returns:
        Path to output file if successful, None otherwise
    """
    print(f"Iniciando alineamiento {algorithm_name} para {os.path.basename(input_fasta)}")
    if os.path.exists(input_fasta) and os.path.getsize(input_fasta) > 0:
        # Get command for the specific alignment algorithm
        cmd = get_command_fn(binary, input_fasta, output_fasta)
        # Handle special case for MAFFT which writes to stdout
        if algorithm_name == "MAFFT":
            with open(output_fasta, "w") as outfile:
                result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
        else:
            result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f"Error ejecutando {algorithm_name} para {input_fasta}:")
            logging.error(result.stderr)
            return None

        print(f"{algorithm_name} ejecutado correctamente, resultados en: {output_fasta}")
        return output_fasta
    else:
        logging.warning(f"No se encontró el archivo {input_fasta}, salteando el alineaminto {algorithm_name}")
        return None

# Command generator functions for each alignment tool
def get_muscle_command(binary, input_fasta, output_fasta):
    cmd = [binary]
    if "5" in binary or "-align" in subprocess.run([binary, "-h"],
                                                    capture_output=True, text=True).stdout:
        cmd.extend(["-align", input_fasta, "-output", output_fasta])
    else:  # Older muscle versions
        cmd.extend(["-in", input_fasta, "-out", output_fasta])
    return cmd

def get_mafft_command(binary, input_fasta, output_fasta):
    return [binary, "--auto", input_fasta]

def get_clustalo_command(binary, input_fasta, output_fasta):
    return [binary, "-i", input_fasta, "-o", output_fasta, "--force", "--outfmt=fasta"]

def run_muscle(input_fasta, output_fasta):
    return run_alignment("MUSCLE", muscle_bin, input_fasta, output_fasta, get_muscle_command)

def run_mafft(input_fasta, output_fasta):
    return run_alignment("MAFFT", mafft_bin, input_fasta, output_fasta, get_mafft_command)

def run_clustalo(input_fasta, output_fasta):
    return run_alignment("Clustal Omega", clustalo_bin, input_fasta, output_fasta, get_clustalo_command)


if __name__ == "__main__":
    # 1. Define los directorios de entrada y salida
    current_dir = os.path.dirname(os.path.abspath(__file__))
    fastas_folder = os.path.join(current_dir, "..", "output", "fastas")
    output_muscle = os.path.join(current_dir, "output_muscle")
    output_mafft = os.path.join(current_dir, "output_mafft")
    output_clustalo = os.path.join(current_dir, "output_clustalo")

    print("Creando directorios de salida...")
    os.makedirs(output_muscle, exist_ok=True)
    os.makedirs(output_mafft, exist_ok=True)
    os.makedirs(output_clustalo, exist_ok=True)
    print("Directorios creados.")

    # Archivos temporales para los alineamientos
    input_dna = "../fastas/lep-sequence_protein_sequences.fasta"
    input_prot = "../fastas/lep-sequence_protein_sequences.fasta"

    # 2. Ejecuta MUSCLE, MAFFT y CLUSTALO para proteínas y ADN
    out_prot_muscle = os.path.join(output_muscle, "aln_proteins.fasta")
    out_dna_muscle = os.path.join(output_muscle, "aln_dna.fasta")
    out_prot_mafft = os.path.join(output_mafft, "aln_proteins.fasta")
    out_dna_mafft = os.path.join(output_mafft, "aln_dna.fasta")
    out_prot_clustalo = os.path.join(output_clustalo, "aln_proteins.fasta")
    out_dna_clustalo = os.path.join(output_clustalo, "aln_dna.fasta")

    prot_aln_muscle = run_muscle(input_prot, out_prot_muscle)
    dna_aln_muscle = run_muscle(input_dna, out_dna_muscle)
    prot_aln_mafft = run_mafft(input_prot, out_prot_mafft)
    dna_aln_mafft = run_mafft(input_dna, out_dna_mafft)
    prot_aln_clustalo = run_clustalo(input_prot, out_prot_clustalo)
    dna_aln_clustalo = run_clustalo(input_dna, out_dna_clustalo)

    # 3. Imprime los alineamientos
    for metodo, prot_aln, dna_aln in [
        ("MUSCLE", prot_aln_muscle, dna_aln_muscle),
        ("MAFFT", prot_aln_mafft, dna_aln_mafft),
        ("CLUSTALO", prot_aln_clustalo, dna_aln_clustalo)
    ]:
        for label, aln_file in [("PROTEÍNAS", prot_aln), ("ADN", dna_aln)]:
            if aln_file and os.path.exists(aln_file):
                print(f"\nAlineamiento {metodo} ({label}):")
                alignment = AlignIO.read(aln_file, "fasta")
                print(f"Número de secuencias alineadas: {len(alignment)}")
                print(f"Longitud del alineamiento: {alignment.get_alignment_length()}")
                for record in alignment:
                    print(f">{record.id}\n{record.seq}\n")
