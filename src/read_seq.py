from Bio import GenBank, SeqIO
import numpy as np


def main():
    with open("../data/lep-sequence.gb") as handle:
        record = SeqIO.read(handle, format="gb")
        print(record.seq)


if __name__ == "__main__":
    main()
