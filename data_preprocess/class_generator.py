'''
Module to generate class labels and balanced multi-FASTA files
for PTM and non-PTM datasets.

The module downsamples the majority class to match a predefined
imbalance factor, then saves the sequences and corresponding
class labels.
'''

import random
import numpy as np
from Bio import SeqIO


def main(positive, negative, output, db='', factor=1.0):
    '''
    Generate balanced multi-FASTA and class label files from
    positive and negative multi-FASTA datasets.

    Args:
        positive (str): Path to the positive multi-FASTA file.
        negative (str): Path to the negative multi-FASTA file.
        output (str): Directory where output files will be saved.
        db (str): Optional prefix for output files.
        factor (float): Ratio of majority to minority class (default 1.0).
    '''
    lens_pos = [len(record.seq) for record in SeqIO.parse(positive, 'fasta')]
    cutoff = np.percentile(lens_pos,95)
    records_pos = [str(record.seq) for record in SeqIO.parse(positive, 'fasta') if len(record.seq) <= cutoff]
    print(f'{positive} len after cutoff {len(records_pos)}')
    records_neg = [str(record.seq) for record in SeqIO.parse(negative, 'fasta') if len(record.seq) <= cutoff]
    print(f'{negative} len after cutoff {len(records_neg)}')

    if len(records_pos) < len(records_neg):
        random.shuffle(records_neg)
        records_neg = records_neg[:int(len(records_pos) * factor)]
    else:
        random.shuffle(records_pos)
        records_pos = records_pos[:int(len(records_neg) * factor)]

    with open(f'{output}/{db}classes.txt', 'w') as file:
        for _ in records_pos:
            file.write('1\n')
        for _ in records_neg:
            file.write('0\n')

    with open(f'{output}/{db}seqs.fasta', 'w') as file:
        i = 1
        for record in records_pos:
            file.write(f'>Seq{i}\n{record}\n')
            i += 1
        for record in records_neg:
            file.write(f'>Seq{i}\n{record}\n')
            i += 1
