# merge existing FASTA files (from databases) and remove duplicate sequences

import os
from Bio import SeqIO
import disease
import matplotlib.pyplot as plt


def main(dir, output):
    merged_file = os.path.join(output, 'merged.fasta')
    if os.path.exists(merged_file):
        os.remove(merged_file)

    merged3000_file = os.path.join(output, 'merged3000.fasta')
    if os.path.exists(merged_file):
        os.remove(merged_file)

    records = []
    new_records = []
    rec_without_disease = []
    labels = []
    all_lens = []

    for filename in os.listdir(dir):
        filepath = os.path.join(dir, filename)
        with_d, without_d = disease.main(filepath)
        rec_without_disease.append(without_d)
        for record in with_d:
            new_records.append(record)
        all_records = records + new_records
        all_records = list(set(all_records))
        # removed = (len(records) + len(new_records)) - len(all_records)
        records = all_records
        # print(f'{filename}: {len(new_records)} sequences. {removed} were already there')
        new_records = []
        labels.append(filename.split('.')[0])
        lens = [len(record.seq) for record in SeqIO.parse(filepath, 'fasta')]
        all_lens.append(lens)
        

    removed_duplicates = []
    seqs = []
    for rec in records:
        if rec[3] not in seqs:
            seqs.append(rec[3])
            removed_duplicates.append(rec)
        else:
            for i, sublist in enumerate(removed_duplicates):
                if rec[3] == sublist[3]:
                    mims = sublist[2] + rec[2]
                    removed_duplicates.pop(i)
                    removed_duplicates.append([rec[0], rec[1], mims, rec[3]])

    lens = []

    with open(merged_file, 'w') as merged:
        for rec in removed_duplicates:
            merged.write(f'>{rec[0]}|{rec[1]}{rec[2]}\n{rec[3]}\n')
            lens.append(len(rec[3]))
        for recs in rec_without_disease:
            for r in recs:
                if r.seq not in seqs:
                    merged.write(f'>{r.description}\n{r.seq}\n')
                    seqs.append(r.seq)
                    lens.append(len(r.seq))
    
    with open(merged3000_file, 'w') as merged:
        for rec in removed_duplicates:
            if(len(rec[3]) <= 3000):
                merged.write(f'>{rec[0]}|{rec[1]}{rec[2]}\n{rec[3]}\n')
        for recs in rec_without_disease:
            for r in recs:
                if r.seq not in seqs and len(r.seq) <= 3000:
                        merged.write(f'>{r.description}\n{r.seq}\n')
                        seqs.append(r.seq)

    fig, axes = plt.subplots(1, 3)

    axes[0].hist(lens, bins=100)
    axes[0].set_title('All sequences')

    axes[1].hist(lens, bins=100, range=[0, 5000])
    axes[1].set_title('Length <= 5000')

    axes[2].hist(lens, bins=100, range=[0, 3000])
    axes[2].set_title('Length <= 3000')

    plt.tight_layout()
    plt.savefig(f'{output}/seq_lens_histo.jpg')
    plt.clf()

    fig, ax = plt.subplots()
    ax.set_ylabel('Sequence length')

    bplot = ax.boxplot(all_lens, tick_labels=labels)
    plt.tight_layout()
    plt.savefig(f'{output}/seq_lens_boxp.jpg')
    plt.clf()

    print(f'{output} Merged succesfully\n')
    disease.vis_disease(output)


if __name__ == '__main__':
    main()
