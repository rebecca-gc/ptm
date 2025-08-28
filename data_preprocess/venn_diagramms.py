'''
Provides visualization of PTM datasets and their disease associations using Venn diagrams.

Functions:
- venn_ptm: Visualizes overlap of sequences with different PTMs.
- venn_disease_ptm: Visualizes overlap of disease-associated sequences across PTMs.
- venn_mim_ids: Visualizes overlap of sequences associated with specific MIM IDs.
- main: Runs all Venn diagram visualizations for a given PTM directory.
'''

import os
import re
import matplotlib.pyplot as plt
from Bio import SeqIO
from venn import venn


def venn_ptm(ptms_dir):
    '''
    Generates a Venn diagram showing the overlap of sequences that have different PTMs.

    Args:
        ptms_dir (str): Path to the directory containing PTM subdirectories with merged multi-FASTA files.
    
    Saves:
        'data/venn_ptm.jpg' containing the Venn diagram of PTM overlaps.
    '''
    sets = []
    names = []
    for d in os.listdir(ptms_dir):
        if os.path.isdir(os.path.join(ptms_dir, d)):
            names.append(d)
            seq_set = set()
            filepath = os.path.join(ptms_dir, d, 'merged.fasta')
            for record in SeqIO.parse(filepath, 'fasta'):
                seq_set.add(record.seq)
            sets.append(seq_set)
    dataset = dict(zip(names, sets))
    venn(dataset)
    plt.savefig('data/venn_ptm.jpg')
    plt.clf()


def venn_disease_ptm(ptms_dir):
    '''
    Generates a Venn diagram showing the overlap of sequences associated with diseases across PTMs.

    Args:
        ptms_dir (str): Path to the directory containing PTM subdirectories with merged multi-FASTA files.
    
    Saves:
        'data/venn_ptm_disease.jpg' containing the Venn diagram of disease-associated PTM sequences.
    '''
    diseases = []
    names = []
    for d in os.listdir(ptms_dir):
        if os.path.isdir(os.path.join(ptms_dir, d)):
            names.append(d)
            disease_set = set()
            filepath = os.path.join(ptms_dir, d, 'merged.fasta')
            for record in SeqIO.parse(filepath, 'fasta'):
                if 'MIM:' in record.description:
                    disease_set.add(record.seq)
            diseases.append(disease_set)
    dataset = dict(zip(names, diseases))
    venn(dataset)
    plt.savefig('data/venn_ptm_disease.jpg')
    plt.clf()


def venn_mim_ids(ptms_dir):
    '''
    Generates Venn diagrams showing the overlap of sequences associated with specific MIM IDs.
    Also generates a version excluding Parkinson's disease for clarity.

    Args:
        ptms_dir (str): Path to the directory containing PTM subdirectories with merged multi-FASTA files.

    Saves:
        'data/venn_disease_overlap5.jpg': all diseases
        'data/venn_disease_overlap4.jpg': all diseases excluding Parkinson's
    '''
    omim_dir = 'local_data/omim'
    seqs = []
    names = []

    for filename in os.listdir(omim_dir):
        if filename.endswith('tsv'):
            mim_ids = []
            seq_set = set()
            names.append(filename.split('.')[0])
            filepath = os.path.join(omim_dir, filename)

            with open(filepath) as file:
                for line in file:
                    if line[0] in '#*%':
                        mim_ids.append(line.split()[0][1:])

            for d in os.listdir(ptms_dir):
                if os.path.isdir(os.path.join(ptms_dir, d)):
                    filepath = os.path.join(ptms_dir, d, 'merged.fasta')
                    for record in SeqIO.parse(filepath, 'fasta'):
                        diseases_in_record = re.findall(r'MIM:(\d+)', record.description)
                        for dis in diseases_in_record:
                            if dis in mim_ids:
                                seq_set.add(record.seq)
            seqs.append(seq_set)

    dataset = dict(zip(names, seqs))
    venn(dataset)
    plt.savefig('data/venn_disease_overlap5.jpg')
    plt.clf()

    dataset = {name: seq for name, seq in zip(names, seqs) if name != 'parkinson'}
    venn(dataset)
    plt.savefig('data/venn_disease_overlap4.jpg')
    plt.clf()


def main(ptms_dir):
    '''
    Executes all Venn diagram visualizations for a given PTM directory.

    Args:
        ptms_dir (str): Path to the directory containing PTM subdirectories with merged multi-FASTA files.
    '''
    venn_ptm(ptms_dir)
    venn_disease_ptm(ptms_dir)
    venn_mim_ids(ptms_dir)
