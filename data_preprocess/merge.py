# merge existing FASTA files (from databases) and remove duplicate sequences

import os
from Bio import SeqIO


GLYCO_DIR = '../data/glycosylation'
S_NITRO_DIR = '../data/s_nitrosylation'
ACET_DIR = '../data/acetylation'


def main():
    dirs = [GLYCO_DIR, S_NITRO_DIR, ACET_DIR]
    for d in dirs:
        merged_file = os.path.join(d, "merged.fasta")
        if os.path.exists(merged_file):
            os.remove(merged_file)

        seqs = []
        new_seqs = []

        for filename in os.listdir(d):
            if filename.endswith(".fasta") and filename != "merged.fasta" and filename != "seqs.fasta":
                filepath = os.path.join(d, filename)
                for record in SeqIO.parse(filepath, "fasta"):
                    new_seqs.append(record.seq)
                all_seqs = seqs + new_seqs
                all_seqs = list(set(all_seqs))
                removed = (len(seqs) + len(new_seqs)) - len(all_seqs)
                seqs = all_seqs
                print(f"{filename}: {len(new_seqs)} sequences. {removed} were already there")
                new_seqs = []

        with open(merged_file, "w") as merged:
            i = 1
            for seq in seqs:
                merged.write(f">Seq{i}\n{seq}\n")
                i += 1

        print("Merged succesfully\n")


if __name__ == '__main__':
    main()
