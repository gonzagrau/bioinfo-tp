from Bio import GenBank, SeqIO
from Bio.Seq import translate
import numpy as np
import re

cds_pattern = r"ATG(?:...)*?(?:ATT|ATC|ACT)"


def main():
    with open("../data/lep-sequence.gb") as handle:
        record = SeqIO.read(handle, format="gb")
        adn_seq = record.seq

    frames = [adn_seq[i:] for i in range(3)] + [
        adn_seq.reverse_complement()[i:] for i in range(3)
    ]

    for i, frame in enumerate(frames):
        print(f"Frame {i+1}: {frame}\n\n")
        cds_list = re.findall(cds_pattern, str(frame))
        for j, cds_seq in enumerate(cds_list):
            translation = translate(cds_seq)
            print(f"Translation {i+1, j}: {translation}\n\n")


if __name__ == "__main__":
    main()
