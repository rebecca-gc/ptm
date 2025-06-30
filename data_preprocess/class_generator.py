import random
from Bio import SeqIO


NO_PTM_FILES = ['data/no_ptm/filtered_no_glyco.fasta',
                'data/no_ptm/filtered_no_s_nitro.fasta',
                'data/no_ptm/filtered_no_acet.fasta',
                'data/no_ptm/filtered_no_methyl.fasta']

DIRS = ['data/glycosylation',
        'data/s_nitrosylation',
        'data/acetylation',
        'data/methylation']


def generator(positive, negative, factor=1.0):
    records_pos = [str(record.seq) for record in SeqIO.parse(f'{positive}/merged.fasta', 'fasta')]
    records_neg = [str(record.seq) for record in SeqIO.parse(negative, 'fasta')]

    if len(records_pos) < len(records_neg):
        random.shuffle(records_neg)
        records_neg = records_neg[:int(len(records_pos)*factor)]
    else:
        random.shuffle(records_pos)
        records_pos = records_pos[:int(len(records_neg)*factor)]

    with open(f'{positive}/classes.txt', 'w') as file:
        for _ in records_pos:
            file.write('1\n')
        for _ in records_neg:
            file.write('0\n')

    print('\nSuccessfully saved classes.txt')

    with open(f'{positive}/seqs.fasta', 'w') as file:
        i = 1
        for record in records_pos:
            file.write(f'>Seq{i}\n{record}\n')
            i += 1
        for record in records_neg:
            file.write(f'>Seq{i}\n{record}\n')
            i += 1

    print('Successfully saved seqs.fasta\n')


def main():
    for i, positive in enumerate(DIRS):
        generator(positive,NO_PTM_FILES[i])


if __name__ == '__main__':
    main()
