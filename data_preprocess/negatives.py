"""
Module for downloading and preparing negative (non-PTM) protein sequences.

This script downloads protein sequences from SwissProt that do NOT have
the specified PTMs and filters out sequences that overlap with positive
PTM datasets to avoid false negatives.

Directory structure:
- NO_PTM_DIR: stores the negative sequences
- DIRS: paths to PTM directories for comparison
"""

import os

import requests
from Bio import SeqIO

NO_PTM_DIR = os.path.join("data", "no_ptm")

DIRS = [
    os.path.join("data", "ptms", "glycosylation"),
    os.path.join("data", "ptms", "s_nitrosylation"),
    os.path.join("data", "ptms", "acetylation"),
    os.path.join("data", "ptms", "methylation"),
]


def swiss_prot(files):
    """
    Download negative sequences from SwissProt for each PTM.

    Each URL queries UniProt for human reviewed sequences excluding
    the given PTM.

    Args:
        files (list of str): filenames to save downloaded sequences
    """
    urls = [
        "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+NOT+%28keyword%3AKW-0325%29%29",
        "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+NOT+%28keyword%3AKW-0702%29%29",
        "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+NOT+%28keyword%3AKW-0007%29%29",
        "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+NOT+%28keyword%3AKW-0488%29%29",
    ]

    for i, url in enumerate(urls):
        filepath = os.path.join(NO_PTM_DIR, files[i])

        with requests.get(url, stream=True, timeout=10) as request:
            request.raise_for_status()
            with open(filepath, "wb") as f:
                for chunk in request.iter_content(chunk_size=2**20):
                    f.write(chunk)


def filter_false_negatives(files):
    """
    Remove sequences from negative datasets that also exist in positive PTM datasets.

    This avoids potential false negatives when training models.

    Args:
        files (list of str): filenames of negative sequences to filter
    """
    for i, file in enumerate(files):
        compare_path = os.path.join(DIRS[i], "merged.fasta")
        target_path = os.path.join(NO_PTM_DIR, file)

        compare_seqs = {
            str(record.seq) for record in SeqIO.parse(compare_path, "fasta")
        }
        target_seqs = {
            str(record.seq) for record in SeqIO.parse(target_path, "fasta")
        }

        common_seqs = compare_seqs & target_seqs
        print(DIRS[i])
        print(f"Count of common sequences: {len(common_seqs)}")

        filtered_path = os.path.join(NO_PTM_DIR, f"filtered_{file}")
        if os.path.exists(filtered_path):
            os.remove(filtered_path)

        with open(filtered_path, "w") as filtered:
            for record in SeqIO.parse(target_path, "fasta"):
                if record.seq not in common_seqs:
                    SeqIO.write(record, filtered, "fasta")

        filtered_seqs = {
            str(record.seq) for record in SeqIO.parse(filtered_path, "fasta")
        }
        print(
            f"PTM: {len(compare_seqs)}, NO_PTM: {len(filtered_seqs)}, before Filter: {len(target_seqs)}"
        )
        print(
            f"percent shared/filtered {len(filtered_seqs) / len(target_seqs)} \n"
        )


def main():
    """
    Main function to download negative sequences and filter false negatives.
    """
    os.makedirs(NO_PTM_DIR, exist_ok=True)

    files = [
        "no_glycosylation.fasta",
        "no_s_nitrosylation.fasta",
        "no_acetylation.fasta",
        "no_methylation.fasta",
    ]

    swiss_prot(files)
    print("Downloads of negative sequences successful\n")

    filter_false_negatives(files)
    print("\nFiltered out all false negatives")
