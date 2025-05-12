# Generates classes.txt and seqs.fasta for two input fasta files

import os
import argparse
from Bio import SeqIO

def generator(positive, negative):

    count_pos = sum(1 for _ in SeqIO.parse(positive, "fasta"))
    count_neg = sum(1 for _ in SeqIO.parse(negative, "fasta"))

    count = count_pos if count_pos < count_neg else count_neg

    with open("classes.txt", "w") as file:
        i = 1
        for _ in SeqIO.parse(positive, "fasta"):
            file.write("1\n")
            i += 1
            if i == count:
                break
        for _ in SeqIO.parse(negative, "fasta"):
            file.write("0\n")
            i += 1
            if i == 2*count-1:
                break
    
    with open("classes.txt", "rb+") as file:
        file.seek(-1, 2)
        file.truncate()

    print('\nSuccessfully saved classes.txt')

    with open("seqs.fasta", "w") as file:
        i = 1
        for record in SeqIO.parse(positive, "fasta"):
            file.write(f">Seq{i}\n{record.seq}\n")
            i += 1
            if i == count:
                break
        for record in SeqIO.parse(negative, "fasta"):
            file.write(f">Seq{i}\n{record.seq}\n")
            i += 1
            if i == 2*count-1:
                break
    
    with open("seqs.fasta", "rb+") as file:
        file.seek(-1, 2)
        file.truncate()

    print('\nSuccessfully saved seqs.fasta\n')

def main():
    positive_error = "Input file path for positive class is bad or the file does not exist"
    negative_error = "Input file path for negative class is bad or the file does not exist"

    parser = argparse.ArgumentParser(description="Generates classes.txt and seqs.fasta from two FASTA input files")
    parser.add_argument("--positive", help="Path to FASTA file with positive samples", required=True)
    parser.add_argument("--negative", help="Path to FASTA file with negative samples", required=True)

    args = parser.parse_args()

    if not os.path.exists(args.positive):
        parser.error(positive_error)
    
    if not os.path.exists(args.negative):
        parser.error(negative_error)

    generator(args.positive,args.negative)

if __name__ == '__main__':
    main()