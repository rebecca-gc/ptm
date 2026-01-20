"""
Main pipeline for PTM dataset preparation and preprocessing.

This script:
1. Downloads raw PTM datasets.
2. Merges datasets from multiple databases.
3. Generates Venn diagrams for overlap visualization.
4. Filters and generates negative datasets.
5. Clusters sequences to reduce redundancy.
6. Generates class labels for positives and negatives.
7. Computes and visualizes sequence length distributions (with 95th-percentile cutoff).
8. Prepares directories for downstream feature importances.
"""

import os

import class_generator
import cluster
import disease
import download_all
import matplotlib.pyplot as plt
import merge
import negatives
import numpy as np
import venn_diagrams
from Bio import SeqIO


def main():
    """Run the full PTM preprocessing pipeline."""

    download_all.main()

    ptms_dir = "data/ptms"

    for ptm in os.listdir(ptms_dir):
        dir_path = os.path.join(ptms_dir, ptm)
        if os.path.isdir(dir_path):
            databases_path = os.path.join(dir_path, "databases")
            merge.main(databases_path, dir_path)

    disease.disease_stacked(ptms_dir)
    venn_diagrams.main(ptms_dir)

    negatives.main()

    for ptm in os.listdir(ptms_dir):
        dir_path = os.path.join(ptms_dir, ptm)
        if os.path.isdir(dir_path):
            databases_path = os.path.join(dir_path, "databases")
            labels = []
            all_lens100 = []
            all_lens95 = []
            for db in os.listdir(databases_path):
                if db.endswith(".fasta"):
                    filepath = os.path.join(databases_path, db)
                    labels.append(db.split(".")[0])
                    lens_100 = [
                        len(record.seq)
                        for record in SeqIO.parse(filepath, "fasta")
                    ]
                    all_lens100.append(lens_100)
                    cutoff = np.percentile(lens_100, 95)
                    lens_95 = [y for y in lens_100 if y <= cutoff]
                    all_lens95.append(lens_95)
                    print(
                        f"{filepath} {max(lens_100)} {max(lens_95)} len: {len(lens_100)}"
                    )
                    short_seqs = [
                        str(record.seq)
                        for record in SeqIO.parse(filepath, "fasta")
                        if len(record.seq) <= 10
                    ]
                    print(short_seqs)

            fig, ax = plt.subplots(figsize=(5, 6))
            ax.set_ylabel("Sequence length", fontsize=14)
            bplot = ax.boxplot(all_lens95, tick_labels=labels)
            plt.xticks(rotation=45)
            ax.tick_params(axis="both", labelsize=14)
            plt.tight_layout()
            plt.savefig(f"{dir_path}/seq_lens_boxp_95.pdf")
            plt.clf()

            fig, ax = plt.subplots(figsize=(5, 6))
            ax.set_ylabel("Sequence length", fontsize=14)
            bplot = ax.boxplot(all_lens100, tick_labels=labels)
            plt.xticks(rotation=45)
            ax.tick_params(axis="both", labelsize=14)
            plt.tight_layout()
            plt.savefig(f"{dir_path}/seq_lens_boxp.pdf")
            plt.clf()


            merged = os.path.join(dir_path, 'merged.fasta')
            cluster.main(merged)

            clustered = os.path.join(dir_path, "clustered.fasta")
            class_generator.main(
                clustered,
                f"data/no_ptm/clustered_no_{ptm}.fasta",
                dir_path,
                factor=1,
            )


            unique = sum(1 for _ in SeqIO.parse(merged, "fasta"))

            clustered_count = sum(1 for _ in SeqIO.parse(clustered, "fasta"))
            uni_clus = clustered_count/unique

            clustered_neg = f'data/no_ptm/clustered_no_{ptm}.fasta'
            neg = sum(1 for _ in SeqIO.parse(clustered_neg, "fasta"))

            print(f'{ptm} unique {unique} clustered {clustered_count} per {uni_clus} neg_clus {neg} \n')


    os.makedirs('data/feature_importances', exist_ok=True)


if __name__ == "__main__":
    main()
