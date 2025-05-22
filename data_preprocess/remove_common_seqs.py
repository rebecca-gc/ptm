import os
import argparse
from Bio import SeqIO


def remove_common_seqs(compare, target):
    """
    remove_common_seqs filters the non-overlapping sequences from
    the target with the compare FASTA file and outputs them into a new FASTA file
    Args:
        compare (os.path): The path to the FASTA file with which the comparison is done.
        target (os.path): The path to the FASTA file out of which the common sequences are filtered.
    Returns:
        None: None
    """

    compare_seqs = {str(record.seq) for record in SeqIO.parse(compare, "fasta")}
    target_seqs = {str(record.seq) for record in SeqIO.parse(target, "fasta")}

    common_seqs = compare_seqs & target_seqs

    print(f"Count of common sequences: {len(common_seqs)}")

    with open("filtered.fasta", "w") as filtered:
        for record in SeqIO.parse(target, "fasta"):
            if record.seq not in common_seqs:
                SeqIO.write(record, filtered, "fasta")

    print("Removed common sequences and saved as filtered.fasta")

    return None


def main():
    compare_file_error = "Input file path for compare file is bad or the file does not exist"
    target_file_error = "Input file path for target file is bad or the file does not exist"

    parser = argparse.ArgumentParser(description="Filters common sequences of two files")
    parser.add_argument("--compare_file", help="Path to compare FASTA file", required=True)
    parser.add_argument("--target_file", help="Path to target FASTA file", required=True)

    args = parser.parse_args()

    if not os.path.exists(args.compare_file):
        parser.error(compare_file_error)

    if not os.path.exists(args.target_file):
        parser.error(target_file_error)

    remove_common_seqs(args.compare_file, args.target_file)

    return None


if __name__ == '__main__':
    main()
