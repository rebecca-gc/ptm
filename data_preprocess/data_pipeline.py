import os
import download_all
import merge
import venn_diagramms
import negatives
import class_generator
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np


def main():
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
            class_generator.generator(os.path.join(dir_path, 'merged.fasta'),f'data/no_ptm/filtered_no_{dir}.fasta',dir_path,factor=1)
            databases_path = os.path.join(dir_path, 'databases')
            db_seqs_classes = os.path.join(dir_path, 'db_seqs_classes')
            try:
                os.makedirs(db_seqs_classes)
            except FileExistsError:
                print(f"Directory '{db_seqs_classes}' already exists.")
            for db in os.listdir(databases_path):
                if db.endswith('fasta'):
                    class_generator.generator(os.path.join(databases_path, db),f'data/no_ptm/filtered_no_{dir}.fasta',db_seqs_classes,db=f'{dir}-{db.split(".")[0]}_',factor=1)
        
            labels = []
            all_lens = []
            for db in os.listdir(databases_path):
                if db.endswith('.fasta'):
                    filepath = os.path.join(databases_path, db)
                    labels.append(db.split('.')[0])
                    lens_100 = [len(record.seq) for record in SeqIO.parse(filepath, 'fasta')]
                    x = np.percentile(lens_100,95)
                    lens_95 = [y for y in lens_100 if y <= x]
                    all_lens.append(lens_95)

            fig, ax = plt.subplots()
            ax.set_ylabel('Sequence length')

            bplot = ax.boxplot(all_lens, tick_labels=labels)
            plt.tight_layout()
            plt.savefig(f'{dir_path}/seq_lens_boxp_95.jpg')
            plt.clf()

    try:
        os.makedirs('data/feature_importances')
    except FileExistsError:
        print(f'Directory "data/feature_importances" already exists.')


    print('\nEverything worked! :)\n')


if __name__ == '__main__':
    main()
