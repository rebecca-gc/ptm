'''
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
'''

import os
import download_all
import merge
import venn_diagramms
import negatives
import cluster
import class_generator
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO


def main():
    '''Run the full PTM preprocessing pipeline.'''

    download_all.main()

    ptms_dir = 'data/ptms'
    for dir in os.listdir(ptms_dir):
        dir_path = os.path.join(ptms_dir, dir)
        if os.path.isdir(dir_path):
            databases_path = os.path.join(dir_path, 'databases')
            merge.main(databases_path, dir_path)

    venn_diagramms.main(ptms_dir)

    negatives.main()

    for dir in os.listdir(ptms_dir):
        dir_path = os.path.join(ptms_dir, dir)
        if os.path.isdir(dir_path):
            merged = os.path.join(dir_path, 'merged.fasta')
            cluster.main(merged)

            clustered = os.path.join(dir_path, 'clustered.fasta')
            class_generator.generator(clustered, f'data/no_ptm/clustered_no_{dir}.fasta', dir_path, factor=1)

            databases_path = os.path.join(dir_path, 'databases')
            labels = []
            all_lens = []
            for db in os.listdir(databases_path):
                if db.endswith('.fasta'):
                    filepath = os.path.join(databases_path, db)
                    labels.append(db.split('.')[0])
                    lens_100 = [len(record.seq) for record in SeqIO.parse(filepath, 'fasta')]
                    cutoff = np.percentile(lens_100, 95)
                    lens_95 = [y for y in lens_100 if y <= cutoff]
                    all_lens.append(lens_95)

            fig, ax = plt.subplots()
            ax.set_ylabel('Sequence length')
            bplot = ax.boxplot(all_lens, tick_labels=labels)
            plt.tight_layout()
            plt.savefig(f'{dir_path}/seq_lens_boxp_95.jpg')
            plt.clf()

    os.makedirs('data/feature_importances', exist_ok=True)


if __name__ == '__main__':
    main()
