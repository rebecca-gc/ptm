# various options for visualization of the datasets

import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from Bio import SeqIO


GLYCO_DIR = '../data/glycosylation'
S_NITRO_DIR = '../data/s_nitrosylation'
ACET_DIR = '../data/acetylation'
dirs = [GLYCO_DIR, S_NITRO_DIR, ACET_DIR]


def venn_ptm():
    sets = []
    for d in dirs:
        seq_set = set()
        filepath = os.path.join(d, "merged.fasta")
        for record in SeqIO.parse(filepath, "fasta"):
            seq_set.add(record.seq)
        sets.append(seq_set)
    venn3(sets, ('Glycosylation', 'S-Nitrosylation', 'Acetylation'))
    plt.title("Comparative Overlap of Protein Sequences with PTMs")
    plt.savefig("venn_ptm.jpg")
    plt.clf()


def venn_disease_ptm():
    diseases = []
    for d in dirs:
        disease = set()
        no_disease = set()
        all_seqs = set()
        filepath = os.path.join(d, "merged.fasta")
        for record in SeqIO.parse(filepath, "fasta"):
            if "MIM:" in record.description:
                disease.add(record.seq)
            else:
                no_disease.add(record.seq)
            all_seqs.add(record.seq)
        diseases.append(disease)
        venn3([disease,no_disease,all_seqs], ('Disease', 'No Disease', f'{d.split("/")[-1]}'))
        plt.title(f'Comparative Overlap of Protein Sequences with {d.split("/")[-1]}')
        plt.savefig(f'{d}/venn_disease.jpg')
        plt.clf()
    venn3(diseases, ('Glycosylation', 'S-Nitrosylation', 'Acetylation'))
    plt.title(f'Comparative Overlap of Protein Sequences associated with Diseases')
    plt.savefig('venn_disease.jpg')
    plt.clf()


def venn_disease():
    return None


def main():
    venn_ptm()
    venn_disease_ptm()


if __name__ == '__main__':
    main()
