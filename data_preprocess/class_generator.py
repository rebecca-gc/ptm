import random
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt


def generator(positive, negative, output, db='', factor=1.0):
    if '3000' in positive:
        records_pos = [str(record.seq) for record in SeqIO.parse(positive, 'fasta')]
        records_neg = [str(record.seq) for record in SeqIO.parse(negative, 'fasta')]
    else:
        # records_pos = [str(record.seq) for record in SeqIO.parse(positive, 'fasta') if len(record.seq) <= 3000]
        lens_pos = [len(record.seq) for record in SeqIO.parse(positive, 'fasta')]
        x = np.percentile(lens_pos,90)
        records_pos = [str(record.seq) for record in SeqIO.parse(positive, 'fasta') if len(record.seq) <= x]
        records_neg = [str(record.seq) for record in SeqIO.parse(negative, 'fasta') if len(record.seq) <= x]

        fig, ax = plt.subplots()
        ax.set_ylabel('Sequence length')

        lens = [len(record) for record in records_pos]
        bplot = ax.boxplot(lens, tick_labels=[db])
        plt.tight_layout()
        plt.savefig(f'{output}/{db}seq_lens_boxp.jpg')
        plt.clf()    

    if len(records_pos) < len(records_neg):
        random.shuffle(records_neg)
        records_neg = records_neg[:int(len(records_pos)*factor)]
    else:
        random.shuffle(records_pos)
        records_pos = records_pos[:int(len(records_neg)*factor)]

    with open(f'{output}/{db}classes.txt', 'w') as file:
        for _ in records_pos:
            file.write('1\n')
        for _ in records_neg:
            file.write('0\n')

    print('\nSuccessfully saved classes.txt')

    with open(f'{output}/{db}seqs.fasta', 'w') as file:
        i = 1
        for record in records_pos:
            file.write(f'>Seq{i}\n{record}\n')
            i += 1
        for record in records_neg:
            file.write(f'>Seq{i}\n{record}\n')
            i += 1

    print('Successfully saved seqs.fasta\n')


def main(positive, negative, output, db, factor):
    generator(positive,negative,output, db, factor)


if __name__ == '__main__':
    main()
