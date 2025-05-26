# merge existing FASTA files and remove duplicate sequences

import os
from Bio import SeqIO


GLYCO_DIR = '/Users/rebeccagrevens/Documents/ptm/data/glycosylation'
S_NITRO_DIR = '/Users/rebeccagrevens/Documents/ptm/data/s_nitrosylation'

def main():
    dirs = [GLYCO_DIR, S_NITRO_DIR]
    for dir in dirs:
        merged_file = os.path.join(dir, "merged.fasta")
        if os.path.exists(merged_file):
            os.remove(merged_file)
        
        seqs = []
        new_seqs = []

        for filename in os.listdir(dir):
            if filename.endswith(".fasta") or filename.endswith(".fa"):
                filepath = os.path.join(dir, filename)
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

        print("\nMerged succesfully")


if __name__ == '__main__':
    main()
