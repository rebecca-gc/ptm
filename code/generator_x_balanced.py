import os
import argparse
import random
from Bio import SeqIO


def generator(positive, negative, factor=1.0):
    """
    generator controls the class imbalance with a factor and
    generates a text file with a 1 for every sample of the positive class
    and a 0 for every sample of the negative class. Additionally, it generates a
    FASTA file for each of the classes
    Args:
        positive (os.path): The path to the FASTA file with samples of the positive class.
        negative (os.path): The path to the FASTA file with samples of the negative class.
        factor (float, optional): The factor that controls the class imbalance. Defaults to 1.0.
    Returns:
        None: None
    """

    records_pos = [str(record.seq) for record in SeqIO.parse(positive, "fasta")]
    records_neg = [str(record.seq) for record in SeqIO.parse(negative, "fasta")]

    if len(records_pos) < len(records_neg):
        random.shuffle(records_neg)
        records_neg = records_neg[:int(len(records_pos)*factor)]
    else:
        random.shuffle(records_pos)
        records_pos = records_pos[:int(len(records_neg)*factor)]

    with open("classes.txt", "w") as file:
        for _ in records_pos:
            file.write("1\n")
        for _ in records_neg:
            file.write("0\n")

    with open("classes.txt", "rb+") as file:
        file.seek(-1, 2)
        file.truncate()

    print('\nSuccessfully saved classes.txt')

    with open("seqs.fasta", "w") as file:
        i = 1
        for record in records_pos:
            file.write(f">Seq{i}\n{record}\n")
            i += 1
        for record in records_neg:
            file.write(f">Seq{i}\n{record}\n")
            i += 1

    with open("seqs.fasta", "rb+") as file:
        file.seek(-1, 2)
        file.truncate()

    print('\nSuccessfully saved seqs.fasta\n')

    return None


def main():
    positive_error = "Input file path for positive class is bad or the file does not exist"
    negative_error = "Input file path for negative class is bad or the file does not exist"
    factor_error = "Factor has to be >= 1"

    parser = argparse.ArgumentParser(description="Generates classes.txt and seqs.fasta from two FASTA input files")
    parser.add_argument("--positive", help="Path to FASTA file with positive samples", required=True)
    parser.add_argument("--negative", help="Path to FASTA file with negative samples", required=True)
    parser.add_argument("--factor", help="Factor for how many more samples of the larger class may be present", required=False, default=1.0)

    args = parser.parse_args()

    if not os.path.exists(args.positive):
        parser.error(positive_error)

    if not os.path.exists(args.negative):
        parser.error(negative_error)

    if float(args.factor) < 1.0:
        parser.error(factor_error)

    generator(args.positive, args.negative, float(args.factor))

    return None


if __name__ == '__main__':
    main()
