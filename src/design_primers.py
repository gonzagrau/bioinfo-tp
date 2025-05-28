import sys, json, os
from Bio import SeqIO
import primer3

def load_config(file: str) -> dict:
    """Load configuration from a JSON file."""	
    with open(file) as f:
        return json.load(f)

def load_seq(fastafile: str) -> tuple:
    """Load a sequence from a FASTA file. returns the sequence and its ID."""
    rec = next(SeqIO.parse(fastafile, "fasta"))
    return str(rec.seq), rec.id

def design(cfg,seq: str, id) -> list:
    """Design primers for a given sequence based on configuration.
    it returns the forward and reverse primers, their melting temperatures, and the product size."""
    ## Prepare the arguments for primer3
    seq_args = {
        "SEQUENCE_ID": id,
        "SEQUENCE_TEMPLATE": seq,
    }
    global_args = {
        "PRIMER_TASK":       "generic",
        "PRIMER_NUM_RETURN": cfg["PRIMER_NUM_RETURN"],
        "PRIMER_MIN_SIZE":   cfg["PRIMER_MIN_SIZE"],
        "PRIMER_OPT_SIZE":   (cfg["PRIMER_MIN_SIZE"] + cfg["PRIMER_MAX_SIZE"]) // 2,
        "PRIMER_MAX_SIZE":   cfg["PRIMER_MAX_SIZE"],
        "PRIMER_MIN_GC":     cfg["PRIMER_MIN_GC_PERCENT"],
        "PRIMER_MAX_GC":     cfg["PRIMER_MAX_GC_PERCENT"],
        "PRIMER_MAX_TM":     cfg["PRIMER_MAX_TM"],
    }
    if cfg.get("AVOID_GC_AT_3_PRIME"):
        global_args["PRIMER_3_END_GC_CLAMP"] = 0

    res = primer3.bindings.design_primers(
        seq_args=seq_args,
        global_args=global_args
    )

    primers = []
    for i in range(cfg["PRIMER_NUM_RETURN"]):
        primers.append({
            "id": id,
            "forward": res[f"PRIMER_LEFT_{i}_SEQUENCE"],
            "tm_fwd":  res[f"PRIMER_LEFT_{i}_TM"],
            "reverse": res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
            "tm_rev":  res[f"PRIMER_RIGHT_{i}_TM"],
            "product": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
        })
    return primers

if __name__=="__main__":
    "para correr el script, usar: python design_primers.py <input.fasta> <params.json> <output.json>"
    "ejemplo: "
#    python .\src\design_primers.py `
#   .\output\fastas\lep-sequence_XP_005250397.1_cds.fasta `
#   .\src\params.json `
#   .\output\primers
    if len(sys.argv)!=4:
        print("Uso: design_primers.py <input.fasta> <params.json> <output.json>")
        sys.exit(1)

    fasta, cfg_fn, out_fn = sys.argv[1], sys.argv[2], sys.argv[3]
    seq, sid = load_seq(fasta)
    cfg = load_config(cfg_fn)
    primers = design( cfg, seq, sid)

    if out_fn:
        os.makedirs(out_fn, exist_ok=True)

    base = os.path.splitext(os.path.basename(fasta))[0]
    out_fn = os.path.join(out_fn, f"resultados_{base}.json")

   
    with open(out_fn, "w") as f:
        json.dump({sid: primers}, f, indent=2)

    print(f"→ {len(primers)} cebadores guardados en {out_fn}")

  
