from Bio import GenBank, SeqIO
import numpy as np
import re

codon_start = "ATG"

def main():
    with open("../data/lep-sequence.gb") as handle:
        record = SeqIO.read(handle, format="gb")
        adn_seq = record.seq
    
    # print(f"{adn_seq}\n\n")

    frames = [adn_seq[i:] for i in range(3)] + [adn_seq.reverse_complement()[i:] for i in range(3)]
 
    for i, frame in enumerate(frames):
        print(f"Frame {i+1}: {frame}\n\n")
        
        match = re.search(codon_start, str(frame))
        
        if not match:
            print(f"No start codon found in frame {i+1}.")
            continue
        
        cds_seg = frame[match.start():]

        translation = cds_seg.translate()
        print(f"Translation {i+1}: {translation}\n\n")
    




if __name__ == "__main__":
    main()
