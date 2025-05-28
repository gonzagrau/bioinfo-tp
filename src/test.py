import sys, json
from Bio import SeqIO
import primer3

def load_config(fn):
    with open(fn) as f:
        return json.load(f)
    

def load_seq(fastafile):
    """Load a sequence from a FASTA file."""
    rec = next(SeqIO.parse(fastafile, "fasta"))
    return str(rec.seq), rec.id


# def design(seq, cfg):
#     p = {
#         "SEQUENCE_ID": "target",
#         "SEQUENCE_TEMPLATE": seq,
#         "PRIMER_TASK": "generic",
#         "PRIMER_NUM_RETURN": cfg["PRIMER_NUM_RETURN"],
#         "PRIMER_MIN_SIZE": cfg["PRIMER_MIN_SIZE"],
#         "PRIMER_OPT_SIZE": (cfg["PRIMER_MIN_SIZE"] + cfg["PRIMER_MAX_SIZE"]) // 2,
#         "PRIMER_MAX_SIZE": cfg["PRIMER_MAX_SIZE"],
#         "PRIMER_MIN_GC": cfg["PRIMER_MIN_GC_PERCENT"],
#         "PRIMER_MAX_GC": cfg["PRIMER_MAX_GC_PERCENT"],
#         "PRIMER_MAX_TM": cfg["PRIMER_MAX_TM"],
#     }
#     if cfg.get("AVOID_GC_AT_3_PRIME"):
#         p["PRIMER_3_END_GC_CLAMP"] = 0

#     res = primer3.bindings.designPrimers({},p)
#     out = []
#     for i in range(cfg["PRIMER_NUM_RETURN"]):
#         out.append({
#             "forward":   res[f"PRIMER_LEFT_{i}_SEQUENCE"],
#             "tm_fwd":    res[f"PRIMER_LEFT_{i}_TM"],
#             "reverse":   res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
#             "tm_rev":    res[f"PRIMER_RIGHT_{i}_TM"],
#             "product":   res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
#         })
#     return out
def design(seq, cfg,id):
    seq_args = {
        "SEQUENCE_ID":id,
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

    # <-- aquí usamos la API moderna:
    res = primer3.bindings.design_primers(
        seq_args=seq_args,
        global_args=global_args
    )

    primers = []
    for i in range(cfg["PRIMER_NUM_RETURN"]):
        primers.append({
            "forward": res[f"PRIMER_LEFT_{i}_SEQUENCE"],
            "tm_fwd":  res[f"PRIMER_LEFT_{i}_TM"],
            "reverse": res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
            "tm_rev":  res[f"PRIMER_RIGHT_{i}_TM"],
            "product": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
        })
    return primers,res

dic = load_config("params.jason")

#print(dic["PRIMER_NUM_RETURN"] )

seq, sid = load_seq("C:/Users/Usuario\Desktop/bioinfo-tp/output/fastas/lep-sequence_NP_000221.1_cds.fasta")
primers,res = design(seq, dic,sid)

#print(primers)
#i want to have a list of the fields in res dic 
fields = [key for key in res.keys()]
print(f"the fields in res are: {fields}")
#print(f"the res is : \n {res}")