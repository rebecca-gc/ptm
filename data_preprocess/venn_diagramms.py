# various options for visualization of the datasets

import os
import re
import matplotlib.pyplot as plt
from Bio import SeqIO
from venn import venn


def venn_ptm(ptms_dir): # this shows the sequences which have x ptms
    sets = []
    names = []
    for d in os.listdir(ptms_dir):
        if os.path.isdir(os.path.join(ptms_dir,d)):
            names.append(d)
            seq_set = set()
            filepath = os.path.join(ptms_dir, d, 'merged.fasta')
            for record in SeqIO.parse(filepath, 'fasta'):
                seq_set.add(record.seq)
            sets.append(seq_set)
    dataset = {name: set_ for name, set_ in zip(names, sets)}
    venn(dataset)
    plt.savefig('data/venn_ptm.jpg')
    plt.clf()


def venn_disease_ptm(ptms_dir): # of the sequences which are associated with diseases, which ones overlap
    diseases = []
    names = []
    for d in os.listdir(ptms_dir):
        if os.path.isdir(os.path.join(ptms_dir,d)):
            names.append(d)
            disease = set()
            # no_disease = set()
            # all_seqs = set()
            filepath = os.path.join(ptms_dir, d, 'merged.fasta')
            for record in SeqIO.parse(filepath, 'fasta'):
                if 'MIM:' in record.description:
                    disease.add(record.seq)
                # else:
                    # no_disease.add(record.seq)
                # all_seqs.add(record.seq)
            diseases.append(disease)
    dataset = {name: dis for name, dis in zip(names, diseases)}
    venn(dataset)
    plt.savefig('data/venn_ptm_disease.jpg')
    plt.clf()


def venn_mim_ids(ptms_dir):
    omim = 'local_data/omim'
    seqs = []
    names = []

    for filename in os.listdir(omim):
        if filename.endswith('tsv'):
            mim_ids = []
            seq_set = set()
            names.append(filename.split('.')[0])

            filepath = os.path.join(omim, filename)
            with open(filepath) as file:
                for line in file:
                    if line[0] == '#' or line[0] == '*' or line[0] == '%':
                        mim_ids.append(line.split()[0][1:])

            for d in os.listdir(ptms_dir):
                if os.path.isdir(os.path.join(ptms_dir,d)):
                    counter = 0
                    filepath = os.path.join(ptms_dir, d, 'merged.fasta')
                    for record in SeqIO.parse(filepath, 'fasta'):
                        diseases = re.findall(r'MIM:(\d+)', record.description)
                        for dis in diseases:
                            if dis in mim_ids:
                                counter += 1
                                seq_set.add(record.seq)
            seqs.append(seq_set)

    dataset = {name: seq for name, seq in zip(names, seqs)}
    venn(dataset)
    plt.savefig('data/venn_disease_overlap5.jpg')
    plt.clf()

    dataset = {name: seq for name, seq in zip(names, seqs) if name != 'parkinson'}
    venn(dataset)
    plt.savefig('data/venn_disease_overlap4.jpg')
    plt.clf()


def main(ptms_dir):
    venn_ptm(ptms_dir)
    venn_disease_ptm(ptms_dir)
    venn_mim_ids(ptms_dir)


if __name__ == '__main__':
    main()
